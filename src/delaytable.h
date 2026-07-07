/*
 * delaytable.h
 * Reader / evaluator for SFXC binary delay tables (.del files).
 *
 * Parses the format produced by SFXC's generate_delay_model (versions -1, 0
 * and 1, cf. sfxc/src/delay_table_akima.cc) and evaluates the delay with a
 * GSL Akima spline over the 1 s samples, exactly as SFXC does.
 *
 * DelayModel aggregates one DelayTable per station and provides a drop-in
 * replacement for the subset of the mpifxcorr Model API used by datasim
 * (getNumScans, getScanStartSec, getScanEndSec, calculateDelayInterpolator).
 * Delays are converted from the SFXC convention (seconds) to the DiFX
 * convention (microseconds, opposite sign): DELAY(us) = -delay_sec * 1e6.
 */
#ifndef DELAYTABLE_H
#define DELAYTABLE_H

#include <cstddef>
#include <string>
#include <vector>

#include <gsl/gsl_spline.h>

/* One phase centre within a scan; Akima spline over scan-relative time. */
class DelaySource
{
 public:
  DelaySource(const std::string &name, const std::vector<double> &t,
              const std::vector<double> &delay);
  ~DelaySource();

  /* Delay in seconds (SFXC sign) at scan-relative time t (seconds). */
  double delay(double t) const;

  const std::string &name() const { return d_name; }
  double tmin() const { return d_tmin; }
  double tmax() const { return d_tmax; }

 private:
  DelaySource(const DelaySource &);
  DelaySource &operator=(const DelaySource &);

  std::string d_name;
  double d_tmin, d_tmax;
  gsl_interp_accel *d_acc;
  gsl_spline *d_spline;
};

/* One scan: begin/end time plus one DelaySource per phase centre. */
class DelayScan
{
 public:
  std::string name;
  int mjd;              // MJD of the scan
  double scan_start;    // seconds of day of scan begin (padding removed)
  double scan_end;      // seconds of day of scan end (padding removed)
  std::vector<DelaySource *> sources;

  DelayScan() : mjd(0), scan_start(0), scan_end(0) {}
  ~DelayScan();
  double beginAbs() const { return mjd * 86400.0 + scan_start; }
  double endAbs() const { return mjd * 86400.0 + scan_end; }

 private:
  DelayScan(const DelayScan &);
  DelayScan &operator=(const DelayScan &);
};

/* Delay table for a single station, read from an SFXC .del file. */
class DelayTable
{
 public:
  explicit DelayTable(const std::string &path);
  ~DelayTable();

  const std::string &station() const { return d_station; }
  const std::string &path() const { return d_path; }
  size_t numScans() const { return d_scans.size(); }
  const DelayScan &scan(size_t idx) const { return *d_scans[idx]; }

  /* Delay in seconds (SFXC sign) at absolute time abssec = mjd*86400 + sec. */
  double delay(size_t scanidx, size_t srcidx, double abssec) const;

 private:
  DelayTable(const DelayTable &);
  DelayTable &operator=(const DelayTable &);

  std::string d_path;
  std::string d_station;
  int d_version;
  int d_npadding;
  std::vector<DelayScan *> d_scans;
};

/* All stations together; drop-in for the Model subset used by datasim. */
class DelayModel
{
 public:
  /* paths: one SFXC delay table per station, in datastream/antenna order.
   * startmjd/startsec: simulation start epoch (offsettime reference).
   * scan: base scan index; scanindex arguments below are relative to it. */
  DelayModel(const std::vector<std::string> &paths, int startmjd, int startsec,
             int scan = 0);
  ~DelayModel();

  int getNumScans() const;
  /* Scan begin/end in whole seconds relative to (mjd, sec). */
  int getScanStartSec(int scanindex, int mjd, int sec) const;
  int getScanEndSec(int scanindex, int mjd, int sec) const;

  const DelayTable &table(size_t antennaindex) const { return *d_tables[antennaindex]; }
  size_t numStations() const { return d_tables.size(); }

  /* Same semantics as mpifxcorr Model::calculateDelayInterpolator (order 0/1):
   * fills delaycoeffs[1] = delay (us) at offsettime and delaycoeffs[0] =
   * delay rate (us per timespan/numincrements) from a 3-point linear fit.
   * offsettime/timespan in seconds relative to the start epoch. */
  bool calculateDelayInterpolator(int scanindex, double offsettime,
                                  double timespan, int numincrements,
                                  int antennaindex, int scansourceindex,
                                  int order, double *delaycoeffs) const;

 private:
  DelayModel(const DelayModel &);
  DelayModel &operator=(const DelayModel &);

  std::vector<DelayTable *> d_tables;
  double d_startabs;    // startmjd*86400 + startsec
  int d_basescan;       // base scan index (JSON "scan" setting)
};

#endif /* DELAYTABLE_H */
