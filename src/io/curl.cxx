#include <curl/curl.h>
#include <json.hpp>

#include "guidescan.hpp"
#include "io/curl.hpp"

using json = nlohmann::json;
using path = std::filesystem::path;


namespace io {
    size_t write_stream(void *ptr, size_t size, size_t nmemb, std::ofstream *stream) {
      stream->write((char*)ptr, size * nmemb);
      return size * nmemb;
    }

    size_t write_string(void *ptr, size_t size, size_t nmemb, std::string *data) {
      data->append((char*)ptr, size * nmemb);
      return size * nmemb;
    }

    int progress_callback(void* clientp, double dltotal, double dlnow, double ultotal, double ulnow)
    {
      if (dltotal <= 0) return 0;
      printf("Downloading %3.0f%\r", dlnow / dltotal * 100);
      fflush(stdout);
      return 0;
    }

    int download_file(std::string url, std::string outfilename) {

      CURL *curl;
      CURLcode res;
      std::ofstream outfile;

      curl = curl_easy_init();
      if (curl) {
        outfile.open(outfilename, std::ios::binary);

        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_stream);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &outfile);
        curl_easy_setopt(curl, CURLOPT_NOPROGRESS, false);
        curl_easy_setopt(curl, CURLOPT_PROGRESSFUNCTION, progress_callback);

        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
        outfile.close();
        return res;
      }
    }

    int download_json(std::string url, json& json_data) {

      CURL *curl;
      CURLcode res;

      std::string json_string;

      curl = curl_easy_init();
      if (curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_string);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &json_string);

        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
      }

      if ((res != CURLE_OK) || (json_string == "")) {
        std::cout << "Server unreachable. You may want to specify the --download-url option." << std::endl;
        return res;
      }

      json_data = json::parse(json_string);
      return 0;
    }
}