namespace phypp {
namespace thread {
    template<typename W, typename T>
    struct worker_with_workspace {
        lock_free_queue<T> input;
        W                  wsp;
        std::atomic<bool>  shutdown;
        std::thread        impl;

        template<typename F, typename ... Args>
        explicit worker_with_workspace(const F& f, const Args&... args) : wsp(args...),
            shutdown(false), impl([this,f]() {

            T t;
            while (!shutdown) {
                while (input.pop(t)) {
                    f(wsp, t);
                }
            }
        }) {}

        ~worker_with_workspace() {
            join();
        }

        void join() {
            if (impl.joinable()) {
                impl.join();
            }
        }

        uint_t workload() const {
            return input.size();
        }
    };

    template<typename T>
    struct worker_no_workspace {
        lock_free_queue<T> input;
        std::atomic<bool>  shutdown;
        std::thread        impl;

        template<typename F>
        explicit worker_no_workspace(const F& f) : shutdown(false),
            impl([this,f]() {

            T t;
            while (!shutdown) {
                while (input.pop(t)) {
                    f(t);
                }
            }
        }) {}

        ~worker_no_workspace() {
            join();
        }

        void join() {
            if (impl.joinable()) {
                impl.join();
            }
        }

        uint_t workload() const {
            return input.size();
        }
    };

    template<typename T, typename W = void>
    struct worker_pool {
        using worker = typename std::conditional<std::is_same<W, void>::value,
            worker_no_workspace<T>, worker_with_workspace<W,T>>::type;

        std::vector<std::unique_ptr<worker>> workers;
        uint_t last_push = 0;

        worker_pool() = default;

        template<typename F, typename ... Args>
        explicit worker_pool(uint_t nthread, const F& f, const Args&... args) {
            start(nthread, f, args...);
        }

        ~worker_pool() {
            join();
        }

        template<typename F, typename ... Args>
        void start(uint_t nthread, const F& f, const Args&... args) {
            workers.clear();
            workers.reserve(nthread);
            for (uint_t i = 0; i < nthread; ++i) {
                workers.emplace_back(new worker(f, args...));
            }

            last_push = workers.size()-1;
        }

        void join() {
            for (uint_t i : range(workers)) {
                workers[i]->shutdown = true;
            }

            for (uint_t i : range(workers)) {
                workers[i]->join();
            }
        }

        void process(T t) {
            ++last_push;
            if (last_push >= workers.size()) {
                last_push = 0;
            }

            workers[last_push]->input.push(std::move(t));
        }

        void process(uint_t i, T t) {
            workers[i]->input.push(std::move(t));
        }

        void consume_all() {
            bool finished = false;
            while (!finished) {
                finished = true;
                for (uint_t i : range(workers)) {
                    if (workers[i]->workload() > 0) {
                        finished = false;
                        break;
                    }
                }

                sleep_for(1e-6);
            }
        }

        uint_t size() const {
            return workers.size();
        }

        uint_t remaining() const {
            uint_t nleft = 0;
            for (uint_t i : range(workers)) {
                nleft += workers[i]->workload();
            }

            return nleft;
        }
    };
}
}
