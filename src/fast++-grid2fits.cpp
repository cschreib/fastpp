#include <phypp.hpp>

const bool debug = false;

struct file_wrapper {
    std::ifstream in;

    explicit file_wrapper(std::string filename) : in(filename, std::ios::binary) {
        in.exceptions(in.failbit);
    }

    template<typename T>
    void read(T& val) {
        in.read(reinterpret_cast<char*>(&val), sizeof(T));
    }

    template<std::size_t D, typename T>
    void read(vec<D,T>& val) {
        uint_t chunk_size = 2e5;
        uint_t nchunk = ceil(val.size()/float(chunk_size));
        if (nchunk > 3) {
            auto pg = progress_start(val.size());
            uint_t nread = 0;
            for (uint_t i : range(nchunk)) {
                uint_t tsize = (i == nchunk-1 ? val.size() - i*chunk_size : chunk_size);
                in.read(reinterpret_cast<char*>(val.data.data()+i*chunk_size), sizeof(T)*tsize);
                if (!in) {
                    throw std::exception();
                }

                nread += tsize;
                print_progress(pg, nread);
            }
        } else {
            in.read(reinterpret_cast<char*>(val.data.data()), sizeof(T)*val.size());
            if (!in) {
                throw std::exception();
            }
        }
    }

    void read(std::string& t) {
        std::uint8_t n = 0;
        read(n);
        t.resize(n);
        in.read(&t[0], t.size());
    }

    template<typename T>
    void seekg(T i) {
        in.seekg(i);
    }
};

int phypp_main(int argc, char* argv[]) {
    if (argc <= 1) {
        print("usage: fast++-grid2fits file.grid");
        return 0;
    }

    std::string filename = argv[1];
    std::string outfile = file::get_directory(filename)+file::get_basename(filename)+".fits";

    if (debug) {
        note("fast++-grid2fits debug mode");
        note("reading ", filename);
        note("writing ", outfile);
    }

    std::uint32_t ngal = 0, ngrid = 0, nprop = 0;
    vec<1,vec1f> grid;
    vec1u grid_dims;
    vec1s grid_names;
    vec1s prop_names;
    uint_t nmodel = 0;

    std::string main_state = "";
    std::string state = "";
    std::string reason = "";

    if (!file::exists(filename)) {
        error("file ", filename, " doest no exist");
        return 1;
    }

    file_wrapper in(filename);

    if (!in.in.is_open()) {
        error("could not open ", filename);
        return 1;
    }

    if (in.in.eof()) {
        error("file ", filename, " is empty");
        return 1;
    }

    if (debug) {
        note("input file successfully open");
    }

    // uint32: size of header in bytes (unused here)
    std::uint32_t hpos = 0; in.read(hpos);
    // unsigned char: file type (C: chi2 grid, B: best chi2)
    unsigned char ftype;
    in.read(ftype);

    if (ftype == 'C') {
        if (debug) {
            note("this is a chi2 grid binary file");
        }

        vec2f chi2;
        vec3f props;

        try {
            // Read header
            main_state = "read header from file";
            if (debug) note(main_state);

            // Format:
            // uint32: number of galaxies
            state = "reading number of galaxies";
            if (debug) note(state);
            in.read(ngal);

            // Check number
            if (ngal > 1e8) {
                reason = "number of galaxies is > 1e8 ("+to_string(ngrid)+")";
                throw std::exception();
            }

            // uint32: number of properties
            state = "reading number of properties";
            if (debug) note(state);
            in.read(nprop);

            // Check number
            if (nprop > 1000) {
                reason = "number of properties is > 1000 ("+to_string(nprop)+")";
                throw std::exception();
            }

            prop_names.resize(nprop);
            for (uint_t i : range(nprop)) {
                in.read(prop_names[i]);
            }

            // uint32: number of grid axis
            state = "reading number of grid parameters";
            if (debug) note(state);
            in.read(ngrid);

            // Check number
            if (nprop > 1000) {
                reason = "number of grid parameter is > 1000 ("+to_string(ngrid)+")";
                throw std::exception();
            }

            // for each grid axis:
            //     uint32: number of values
            //     float[*]: grid values
            nmodel = 1;
            grid.resize(ngrid);
            grid_dims.resize(ngrid);
            grid_names.resize(ngrid);
            for (uint_t i : range(grid)) {
                state = "reading grid for parameter "+to_string(i);
                if (debug) note(state);
                in.read(grid_names[i]);

                std::uint32_t gsize;
                in.read(gsize);

                grid_dims[i] = gsize;
                nmodel *= grid_dims[i];
                grid[i].resize(gsize);
                in.read(grid[i]);
            }

            print("found ", ngal, " galaxies, with ", nprop, " properties, ", ngrid,
                " grid parameters and ", nmodel, " models");
            print("grid names: ", collapse(grid_names, ", "));
            print("prop names: ", collapse(prop_names, ", "));

            main_state = "read data";

            state = "properties";
            props.resize(nmodel, ngal, 1+nprop);
            in.read(props);

            prepend(prop_names, vec1s{"chi2"});
        } catch (...) {
            error("\n\ncould not ", main_state, " (", state, ")");
            return 1;
        }

        vec3f vgrid(nmodel, ngal, ngrid);
        vec1u idm(ngrid);
        for (uint_t im = 0; im < nmodel; ++im) {
            for (uint_t ig : range(ngrid)) {
                vgrid.safe(im,_,ig) = grid.safe[ig].safe[idm.safe[ig]];
            }

            increment_index_list(idm, grid_dims);
        }

        print("writing ", outfile);

        fits::output_table otbl(outfile);
        otbl.reach_hdu(1);

        for (uint_t i : range(prop_names)) {
            otbl.write_column(prop_names[i], props(_,_,i));
        }
        for (uint_t i : range(grid_names)) {
            otbl.write_column(grid_names[i], vgrid(_,_,i));
        }
    } else if (ftype == 'B') {
        if (debug) {
            note("this is a best chi2 binary file");
        }

        vec1u id_model;
        vec1f chi2;
        vec2f props;

        try {
            // Read header
            main_state = "read header from file";

            // Format:
            // uint32: number of properties
            state = "reading number of properties";
            if (debug) note(state);
            in.read(nprop);

            // Check number
            if (nprop > 1000) {
                reason = "number of properties is > 1000 ("+to_string(nprop)+")";
                throw std::exception();
            }

            prop_names.resize(nprop);
            for (uint_t i : range(nprop)) {
                in.read(prop_names[i]);
            }

            // uint32: number of grid axis
            state = "reading number of grid parameters";
            if (debug) note(state);
            in.read(ngrid);

            // Check number
            if (nprop > 1000) {
                reason = "number of grid parameter is > 1000 ("+to_string(ngrid)+")";
                throw std::exception();
            }

            // for each grid axis:
            //     uint32: number of values
            //     float[*]: grid values
            grid.resize(ngrid);
            grid_dims.resize(ngrid);
            grid_names.resize(ngrid);
            for (uint_t i : range(grid)) {
                state = "reading grid for parameter "+to_string(i);
                if (debug) note(state);
                in.read(grid_names[i]);
                std::uint32_t gsize; in.read(gsize);

                // Check number
                if (gsize > 100000) {
                    reason = "size of grid for "+grid_names[i]+" is > 100000 ("+to_string(gsize)+")";
                    throw std::exception();
                }

                grid_dims[i] = gsize;
                grid[i].resize(gsize);
                in.read(grid[i]);
            }

            main_state = "read data";
            state = "";

            // As a failsafe, compute maximum possible iterations
            uint_t max_iter; {
                auto opos = in.in.tellg();
                in.in.seekg(0, std::ios::end);
                uint_t flen = in.in.tellg() - opos;
                in.in.seekg(opos);
                uint_t data_size = sizeof(uint32_t) + sizeof(float) + sizeof(float)*nprop;
                max_iter = flen/data_size + 2;
            }

            // Try to read as much as we can, writing may have been interrupted in the middle
            // of the computations, and we want to salvage whatever is available
            uint_t iter = 0;
            while (iter < max_iter) {
                uint32_t id;
                float tchi2;
                vec1f p(nprop);

                try {
                    in.read(id);
                    in.read(tchi2);
                    in.read(p);

                    id_model.push_back(id);
                    chi2.push_back(tchi2);
                    props.push_back(p);
                } catch (...) {
                    break;
                }

                ++iter;
            }

            if (iter == max_iter) {
                reason = "maximum number of iteration reached, that's a bug...";
                throw std::exception();
            }

            nmodel = id_model.size();

            print("found ", nprop, " properties, ", ngrid, " grid parameters and ", nmodel, " models");
            print("grid names: ", collapse(grid_names, ", "));
            print("prop names: ", collapse(prop_names, ", "));
        } catch (...) {
            print("\n\n");
            error("could not ", main_state, " (", state, ")");
            if (!reason.empty()) {
                error(reason);
            }
            return 1;
        }

        vec2f vgrid(nmodel, ngrid);

        for (uint_t im : range(nmodel)) {
            vec1u idm(ngrid); {
                uint_t i = id_model[im];
                for (uint_t j : range(ngrid)) {
                    idm[ngrid-1-j] = i % grid_dims[ngrid-1-j];
                    i /= grid_dims[ngrid-1-j];
                }
            }

            for (uint_t ig : range(ngrid)) {
                vgrid.safe(im,ig) = grid.safe[ig].safe[idm.safe[ig]];
            }
        }

        print("writing ", outfile);

        fits::output_table otbl(outfile);
        otbl.reach_hdu(1);

        otbl.write_column("chi2", chi2);
        otbl.write_column("model", id_model);

        for (uint_t i : range(prop_names)) {
            otbl.write_column(prop_names[i], props(_,i));
        }
        for (uint_t i : range(grid_names)) {
            otbl.write_column(grid_names[i], vgrid(_,i));
        }
    } else {
        error("unknown file type '", ftype, "'; are you sure this is a *.grid file from FAST++?");
        return 1;
    }

    return 0;
}
