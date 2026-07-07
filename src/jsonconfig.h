/*
 * jsonconfig.h
 * JSON-based configuration for datasim, replacing the DiFX .input file
 * (mpifxcorr Configuration class).  Getter names mirror the subset of the
 * Configuration API used by datasim; configindex arguments are accepted but
 * ignored.
 *
 * Schema (all frequencies/bandwidths in MHz):
 * {
 *   "start": { "mjd": 60453, "seconds": 43200 },
 *   "duration": 120.0,
 *   "flux_density": 100.0,          // Jy (optional, default 100)
 *   "seed": 48573,                  // optional
 *   "numdivs": 1,                   // optional, time-division parallelisation
 *   "specres_scale": 1,             // optional
 *   "pcal": 0,                      // optional, phase-cal interval
 *   "line_signal": [],              // optional, [freq, width, amplitude]
 *   "scan": 0,                      // optional, delay-table scan index
 *   "stations": [
 *     { "name": "EF",
 *       "sefd": 1000,
 *       "delay_table": "delays/N24L2_Ef.del",
 *       "framebytes": 8032,
 *       "frames_per_sec": 1000,     // optional, checked against derived value
 *       "bands": [ { "freq": 1626.49, "bandwidth": 16.0 } ] }
 *   ]
 * }
 */
#ifndef JSONCONFIG_H
#define JSONCONFIG_H

#include <string>
#include <vector>

class JsonConfig
{
 public:
  explicit JsonConfig(const std::string &path);

  /* --- Configuration-compatible getters (configindex ignored) --- */
  int getStartMJD() const { return d_startmjd; }
  int getStartSeconds() const { return d_startseconds; }
  double getExecuteSeconds() const { return d_duration; }
  int getNumDataStreams() const { return (int)d_stations.size(); }
  std::string getTelescopeName(int antidx) const { return station(antidx).name; }
  int getFrameBytes(int /*configindex*/, int antidx) const { return station(antidx).framebytes; }
  int getFramesPerSecond(int /*configindex*/, int antidx) const { return station(antidx).framespersec; }
  int getDNumRecordedBands(int /*configindex*/, int antidx) const { return (int)station(antidx).freq.size(); }
  double getDRecordedFreq(int /*configindex*/, int antidx, int bandidx) const;
  double getDRecordedBandwidth(int /*configindex*/, int antidx, int bandidx) const;
  int getDLocalRecordedFreqIndex(int /*configindex*/, int /*antidx*/, int bandidx) const { return bandidx; }

  /* --- datasim-specific settings --- */
  const std::string &getDelayTable(int antidx) const { return station(antidx).delaytable; }
  int getSEFD(int antidx) const { return station(antidx).sefd; }
  float getFluxDensity() const { return d_fluxdensity; }
  unsigned int getSeed() const { return d_seed; }
  int getNumDivs() const { return d_numdivs; }
  int getPcal() const { return d_pcal; }
  int getSpecResScale() const { return d_specresscale; }
  int getScan() const { return d_scan; }
  const std::vector<float> &getLineSignal() const { return d_linesignal; }

 private:
  struct Station
  {
    std::string name;
    std::string delaytable;
    int sefd;
    int framebytes;
    int framespersec;
    std::vector<double> freq;        // MHz, per recorded band
    std::vector<double> bandwidth;   // MHz, per recorded band
  };

  const Station &station(int antidx) const;

  int d_startmjd;
  int d_startseconds;
  double d_duration;
  float d_fluxdensity;
  unsigned int d_seed;
  int d_numdivs;
  int d_pcal;
  int d_specresscale;
  int d_scan;
  std::vector<float> d_linesignal;
  std::vector<Station> d_stations;
};

#endif /* JSONCONFIG_H */
