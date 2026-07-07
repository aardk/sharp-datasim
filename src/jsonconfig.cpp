/*
 * jsonconfig.cpp
 * JSON configuration parser for datasim; see jsonconfig.h for the schema.
 */
#include "jsonconfig.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "json.hpp"

using nlohmann::json;

namespace {

const int VDIF_HDR_BYTES = 32;
const int BITS_PER_SAMPLE = 2;

/* Derived VDIF frames per second: Nyquist-sampled real data, 2 bits/sample,
 * all bands of a datastream interleaved in one frame (DiFX convention). */
int deriveFramesPerSec(const std::string &name, double bandwidthMHz,
                       size_t nbands, int framebytes)
{
  double payload = framebytes - VDIF_HDR_BYTES;
  if(payload <= 0)
    throw std::runtime_error("station '" + name + "': framebytes must be > 32");
  double bytespersec =
    nbands * bandwidthMHz * 1e6 * 2.0 * BITS_PER_SAMPLE / 8.0;
  double fps = bytespersec / payload;
  if(fps <= 0.0 || fps != (double)(int)fps)
  {
    std::ostringstream oss;
    oss << "station '" << name << "': frames per second (" << fps
        << ") derived from bandwidth/framebytes is not a positive integer";
    throw std::runtime_error(oss.str());
  }
  return (int)fps;
}

}  // namespace

JsonConfig::JsonConfig(const std::string &path)
  : d_startmjd(0), d_startseconds(0), d_duration(0), d_fluxdensity(100.0f),
    d_seed(48573), d_numdivs(1), d_pcal(0), d_specresscale(1), d_scan(0)
{
  std::ifstream in(path.c_str());
  if(!in)
    throw std::runtime_error("JsonConfig: cannot open '" + path + "'");

  json cfg;
  try
  {
    cfg = json::parse(in, NULL, true, true);   // allow // and /* */ comments
  }
  catch(const json::parse_error &e)
  {
    throw std::runtime_error("JsonConfig: '" + path + "': " + e.what());
  }

  try
  {
    const json &start = cfg.at("start");
    d_startmjd = start.at("mjd").get<int>();
    d_startseconds = start.at("seconds").get<int>();
    d_duration = cfg.at("duration").get<double>();
    if(d_duration <= 0)
      throw std::runtime_error("duration must be > 0");

    d_fluxdensity = cfg.value("flux_density", 100.0);
    d_seed = cfg.value("seed", 48573u);
    d_numdivs = cfg.value("numdivs", 1);
    d_pcal = cfg.value("pcal", 0);
    d_specresscale = cfg.value("specres_scale", 1);
    d_scan = cfg.value("scan", 0);
    if(cfg.contains("line_signal"))
    {
      const json &ls = cfg.at("line_signal");
      if(!ls.is_array() || (ls.size() != 0 && ls.size() != 3))
        throw std::runtime_error("line_signal must be [freq, width, amplitude]");
      for(size_t i = 0; i < ls.size(); i++)
        d_linesignal.push_back(ls[i].get<float>());
    }

    const json &stations = cfg.at("stations");
    if(!stations.is_array() || stations.empty())
      throw std::runtime_error("stations must be a non-empty array");
    for(size_t i = 0; i < stations.size(); i++)
    {
      const json &st = stations[i];
      Station s;
      s.name = st.at("name").get<std::string>();
      s.delaytable = st.at("delay_table").get<std::string>();
      s.sefd = st.at("sefd").get<int>();
      s.framebytes = st.at("framebytes").get<int>();

      const json &bands = st.at("bands");
      if(!bands.is_array() || bands.empty())
        throw std::runtime_error("station '" + s.name +
                                 "': bands must be a non-empty array");
      for(size_t j = 0; j < bands.size(); j++)
      {
        s.freq.push_back(bands[j].at("freq").get<double>());
        s.bandwidth.push_back(bands[j].at("bandwidth").get<double>());
        if(s.bandwidth.back() != s.bandwidth.front())
          throw std::runtime_error("station '" + s.name +
                                   "': all bands must have equal bandwidth");
      }

      s.framespersec = deriveFramesPerSec(s.name, s.bandwidth[0],
                                          s.freq.size(), s.framebytes);
      if(st.contains("frames_per_sec") &&
         st.at("frames_per_sec").get<int>() != s.framespersec)
      {
        std::ostringstream oss;
        oss << "station '" << s.name << "': frames_per_sec ("
            << st.at("frames_per_sec").get<int>() << ") does not match value "
            << "derived from bandwidth/framebytes (" << s.framespersec << ")";
        throw std::runtime_error(oss.str());
      }
      d_stations.push_back(s);
    }
  }
  catch(const json::exception &e)
  {
    throw std::runtime_error("JsonConfig: '" + path + "': " + e.what());
  }
  catch(const std::runtime_error &e)
  {
    throw std::runtime_error("JsonConfig: '" + path + "': " + e.what());
  }
}

const JsonConfig::Station &JsonConfig::station(int antidx) const
{
  if(antidx < 0 || (size_t)antidx >= d_stations.size())
    throw std::out_of_range("JsonConfig: station index out of range");
  return d_stations[antidx];
}

double JsonConfig::getDRecordedFreq(int, int antidx, int bandidx) const
{
  const Station &s = station(antidx);
  if(bandidx < 0 || (size_t)bandidx >= s.freq.size())
    throw std::out_of_range("JsonConfig: band index out of range");
  return s.freq[bandidx];
}

double JsonConfig::getDRecordedBandwidth(int, int antidx, int bandidx) const
{
  const Station &s = station(antidx);
  if(bandidx < 0 || (size_t)bandidx >= s.bandwidth.size())
    throw std::out_of_range("JsonConfig: band index out of range");
  return s.bandwidth[bandidx];
}
