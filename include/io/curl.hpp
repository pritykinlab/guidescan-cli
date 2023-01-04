#include <json.hpp>

using json = nlohmann::json;

namespace io {
    int download_file(std::string url, std::string outfilename);
    json download_json(std::string url);
}
