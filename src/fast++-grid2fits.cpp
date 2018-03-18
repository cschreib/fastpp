#include <phypp.hpp>

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

                nread += tsize;
                print_progress(pg, nread);
            }
        } else {
            in.read(reinterpret_cast<char*>(val.data.data()), sizeof(T)*val.size());
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

    std::uint32_t ngal = 0, ngrid = 0, nprop = 0;
    vec<1,vec1f> grid;
    vec1u grid_dims;
    vec1s grid_names;
    vec1s prop_names;
    uint_t nmodel = 0;

    std::string main_state = "";
    std::string state = "";

    file_wrapper in(filename);

    if (!in.in.is_open()) {
        error("could not open ", filename);
        return 1;
    }
    if (in.in.eof()) {
        error("file ", filename, " is empty");
        return 1;
    }

    // uint32: size of header in bytes (unused here)
    std::uint32_t hpos = 0; in.read(hpos);
    // unsigned char: file type (C: chi2 grid, B: best chi2)
    unsigned char ftype;
    in.read(ftype);

    if (ftype == 'C') {
        vec2f chi2;
        vec3f props;

        try {
            // Read header
            main_state = "read header from file";

            // Format:
            // uint32: number of galaxies
            state = "reading number of galaxies"; in.read(ngal);
            // uint32: number of properties
            state = "reading number of properties"; in.read(nprop);
            prop_names.resize(nprop);
            for (uint_t i : range(nprop)) {
                in.read(prop_names[i]);
            }

            // uint32: number of grid axis
            state = "reading number of grid parameters"; in.read(ngrid);
            // for each grid axis:
            //     uint32: number of values
            //     float[*]: grid values
            nmodel = 1;
            grid.resize(ngrid);
            grid_dims.resize(ngrid);
            grid_names.resize(ngrid);
            for (uint_t i : range(grid)) {
                state = "reading grid for parameter "+to_string(i);
                in.read(grid_names[i]);

                std::uint32_t gsize;
                in.read(gsize);

                grid_dims[i] = gsize;
                nmodel *= grid_dims[i];
                grid[i].resize(gsize);
                in.read(grid[i]);
            }

            print("found ", ngal, " galaxies, with ", nprop, " properties, ", ngrid,
                " parameters and ", nmodel, " models");
            print("grid names: ", collapse(grid_names, ", "));
            print("prop names: ", collapse(prop_names, ", "));

            main_state = "read data";

            state = "properties";
            props.resize(nmodel, ngal, 1+nprop);
            in.read(props);

            prepend(prop_names, vec1s{"chi2"});
        } catch (...) {
            error("\n\ncould not ", main_state, " (", state, ")");
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
        vec1u id_model;
        vec1f chi2;
        vec2f props;

        try {
            // Read header
            main_state = "read header from file";

            // Format:
            // uint32: number of properties
            state = "reading number of properties"; in.read(nprop);
            prop_names.resize(nprop);
            for (uint_t i : range(nprop)) {
                in.read(prop_names[i]);
            }

            // uint32: number of grid axis
            state = "reading number of grid parameters"; in.read(ngrid);
            // for each grid axis:
            //     uint32: number of values
            //     float[*]: grid values
            grid.resize(ngrid);
            grid_dims.resize(ngrid);
            grid_names.resize(ngrid);
            for (uint_t i : range(grid)) {
                state = "reading grid for parameter "+to_string(i);
                in.read(grid_names[i]);
                std::uint32_t gsize; in.read(gsize);
                grid_dims[i] = gsize;
                grid[i].resize(gsize);
                in.read(grid[i]);
            }

            main_state = "read data";

            // Try to read as much as we can, writing may have been interrupted
            while (true) {
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
            }

            nmodel = id_model.size();

            print("found ", nprop, " properties, ", ngrid, " parameters and ", nmodel, " models");
            print("grid names: ", collapse(grid_names, ", "));
            print("prop names: ", collapse(prop_names, ", "));
        } catch (...) {
            error("\n\ncould not ", main_state, " (", state, ")");
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
