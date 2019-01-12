#ifndef FASTPP_IO_HPP
#define FASTPP_IO_HPP

// Helper functions
// ----------------

namespace vif {
namespace file {
    template<typename S, typename T>
    bool write(S& s, const T& t) {
        s.write(reinterpret_cast<const char*>(&t), sizeof(T));
        return !s.fail();
    }

    template<typename R, typename S, typename T>
    bool write_as(S& s, const T& t) {
        R r = t;
        s.write(reinterpret_cast<const char*>(&r), sizeof(R));
        return !s.fail();
    }

    template<typename S>
    bool write(S& s, const std::string& t) {
        write_as<std::uint8_t>(s, t.size());
        s.write(t.c_str(), t.size());
        return !s.fail();
    }

    template<typename S, std::size_t D, typename T>
    bool write(S& s, const vec<D,T>& t) {
        s.write(reinterpret_cast<const char*>(t.data.data()), sizeof(T)*t.size());
        return !s.fail();
    }

    template<typename S, typename T>
    bool read(S& s, T& t) {
        s.read(reinterpret_cast<char*>(&t), sizeof(T));
        return !s.fail();
    }

    template<typename S>
    bool read(S& s, std::string& t) {
        std::uint8_t n = 0;
        if (read(s, n)) {
            t.resize(n);
            s.read(&t[0], t.size());
            return !s.fail();
        } else {
            return false;
        }
    }

    template<typename R, typename S, typename T>
    bool read_as(S& s, T& t) {
        R r;
        if (s.read(reinterpret_cast<char*>(&r), sizeof(R))) {
            t = r;
        }

        return !s.fail();
    }

    template<typename S, std::size_t D, typename T>
    bool read(S& s, vec<D,T>& t) {
        s.read(reinterpret_cast<char*>(t.data.data()), sizeof(T)*t.size());
        return !s.fail();
    }
}
}

#endif
