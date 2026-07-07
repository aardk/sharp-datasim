/*
 * delaytable.cpp
 * Implementation of the SFXC binary delay table reader/evaluator.
 * See delaytable.h and sfxc/src/delay_table_akima.cc for the format.
 */
#include "delaytable.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

const size_t NAME_LEN = 81;      // scan/source name field size
const size_t REC_DOUBLES = 7;    // time, u, v, w, delay, phase, amplitude

std::string decodeName(const char *raw, size_t len)
{
  size_t n = 0;
  while(n < len && raw[n] != '\0') n++;
  std::string s(raw, n);
  size_t b = s.find_first_not_of(" \t");
  size_t e = s.find_last_not_of(" \t");
  if(b == std::string::npos) return std::string();
  return s.substr(b, e - b + 1);
}

}  // namespace

/* ------------------------------------------------------------------ */
/* DelaySource                                                         */
/* ------------------------------------------------------------------ */
DelaySource::DelaySource(const std::string &name, const std::vector<double> &t,
                         const std::vector<double> &delay)
  : d_name(name), d_acc(NULL), d_spline(NULL)
{
  if(t.size() < 5)
    throw std::runtime_error("DelaySource '" + name + "': need at least 5 "
                             "samples for an Akima spline");
  d_tmin = t.front();
  d_tmax = t.back();
  d_acc = gsl_interp_accel_alloc();
  d_spline = gsl_spline_alloc(gsl_interp_akima, t.size());
  gsl_spline_init(d_spline, &t[0], &delay[0], t.size());
}

DelaySource::~DelaySource()
{
  if(d_spline) gsl_spline_free(d_spline);
  if(d_acc) gsl_interp_accel_free(d_acc);
}

double DelaySource::delay(double t) const
{
  if(t < d_tmin || t > d_tmax)
  {
    std::ostringstream oss;
    oss << "DelaySource '" << d_name << "': time " << t
        << " outside spline domain [" << d_tmin << ", " << d_tmax << "]";
    throw std::out_of_range(oss.str());
  }
  return gsl_spline_eval(d_spline, t, d_acc);
}

/* ------------------------------------------------------------------ */
/* DelayScan                                                           */
/* ------------------------------------------------------------------ */
DelayScan::~DelayScan()
{
  for(size_t i = 0; i < sources.size(); i++)
    delete sources[i];
}

/* ------------------------------------------------------------------ */
/* DelayTable                                                          */
/* ------------------------------------------------------------------ */
DelayTable::DelayTable(const std::string &path)
  : d_path(path), d_version(-1), d_npadding(0)
{
  std::ifstream in(path.c_str(), std::ios::binary);
  if(!in)
    throw std::runtime_error("DelayTable: cannot open '" + path + "'");

  int32_t header_size = 0;
  if(!in.read(reinterpret_cast<char *>(&header_size), sizeof(header_size)))
    throw std::runtime_error("DelayTable: '" + path + "' truncated header");
  std::vector<char> header(header_size > 0 ? header_size : 0);
  if(header_size > 0 &&
     !in.read(&header[0], header_size))
    throw std::runtime_error("DelayTable: '" + path + "' truncated header");

  if(header_size >= 4)
    memcpy(&d_version, &header[0], sizeof(int32_t));
  if(d_version < -1 || d_version > 1)
  {
    std::ostringstream oss;
    oss << "DelayTable: '" << path << "' unsupported version " << d_version;
    throw std::runtime_error(oss.str());
  }
  if(d_version >= 1)
  {
    memcpy(&d_npadding, &header[4], sizeof(int32_t));
    d_station = decodeName(&header[8], 3);
  }

  // Read per-source blocks: [scan name (v>=0)] source name, mjd, records
  // terminated by a record with time == 0 && delay == 0.
  char namebuf[NAME_LEN];
  while(in.peek() != EOF)
  {
    std::string scan_name;
    if(d_version >= 0)
    {
      if(!in.read(namebuf, NAME_LEN)) break;
      scan_name = decodeName(namebuf, NAME_LEN);
    }
    if(!in.read(namebuf, NAME_LEN)) break;
    std::string source_name = decodeName(namebuf, NAME_LEN);
    int32_t mjd = 0;
    if(!in.read(reinterpret_cast<char *>(&mjd), sizeof(mjd))) break;

    std::vector<double> times, delays;
    double rec[REC_DOUBLES];
    while(in.read(reinterpret_cast<char *>(rec), sizeof(rec)))
    {
      if(rec[0] == 0.0 && rec[4] == 0.0) break;   // scan terminator
      times.push_back(rec[0]);
      delays.push_back(rec[4]);
    }
    if(times.empty()) continue;

    double scan_start = times.front() + d_npadding;
    double scan_end = times.back() - d_npadding;
    std::vector<double> trel(times.size());
    for(size_t i = 0; i < times.size(); i++)
      trel[i] = times[i] - scan_start;

    // Consecutive blocks sharing the scan begin time are multiple phase
    // centres of the same scan.
    bool same = !d_scans.empty() && d_scans.back()->mjd == mjd &&
                fabs(d_scans.back()->scan_start - scan_start) < 0.5e-3;
    if(!same)
    {
      DelayScan *sc = new DelayScan();
      sc->name = scan_name;
      sc->mjd = mjd;
      sc->scan_start = scan_start;
      sc->scan_end = scan_end;
      d_scans.push_back(sc);
    }
    d_scans.back()->sources.push_back(new DelaySource(source_name, trel, delays));
  }

  if(d_scans.empty())
    throw std::runtime_error("DelayTable: '" + path + "' contains no scans");
}

DelayTable::~DelayTable()
{
  for(size_t i = 0; i < d_scans.size(); i++)
    delete d_scans[i];
}

double DelayTable::delay(size_t scanidx, size_t srcidx, double abssec) const
{
  if(scanidx >= d_scans.size())
    throw std::out_of_range("DelayTable: scan index out of range");
  const DelayScan &sc = *d_scans[scanidx];
  if(srcidx >= sc.sources.size())
    throw std::out_of_range("DelayTable: source index out of range");
  return sc.sources[srcidx]->delay(abssec - sc.beginAbs());
}

/* ------------------------------------------------------------------ */
/* DelayModel                                                          */
/* ------------------------------------------------------------------ */
DelayModel::DelayModel(const std::vector<std::string> &paths, int startmjd,
                       int startsec, int scan)
  : d_startabs(startmjd * 86400.0 + startsec), d_basescan(scan)
{
  if(paths.empty())
    throw std::runtime_error("DelayModel: no delay tables given");
  try
  {
    for(size_t i = 0; i < paths.size(); i++)
      d_tables.push_back(new DelayTable(paths[i]));
  }
  catch(...)
  {
    for(size_t i = 0; i < d_tables.size(); i++)
      delete d_tables[i];
    throw;
  }

  // All stations must agree on the scan layout.
  for(size_t i = 1; i < d_tables.size(); i++)
  {
    if(d_tables[i]->numScans() != d_tables[0]->numScans())
      throw std::runtime_error("DelayModel: '" + paths[i] + "' has a "
                               "different number of scans than '" + paths[0] + "'");
  }
  if(d_basescan < 0 || (size_t)d_basescan >= d_tables[0]->numScans())
    throw std::runtime_error("DelayModel: base scan index out of range");
}

DelayModel::~DelayModel()
{
  for(size_t i = 0; i < d_tables.size(); i++)
    delete d_tables[i];
}

int DelayModel::getNumScans() const
{
  return static_cast<int>(d_tables[0]->numScans());
}

int DelayModel::getScanStartSec(int scanindex, int mjd, int sec) const
{
  const DelayScan &sc = d_tables[0]->scan(d_basescan + scanindex);
  return static_cast<int>(floor(sc.beginAbs() - (mjd * 86400.0 + sec)));
}

int DelayModel::getScanEndSec(int scanindex, int mjd, int sec) const
{
  const DelayScan &sc = d_tables[0]->scan(d_basescan + scanindex);
  return static_cast<int>(floor(sc.endAbs() - (mjd * 86400.0 + sec)));
}

bool DelayModel::calculateDelayInterpolator(int scanindex, double offsettime,
                                            double timespan, int numincrements,
                                            int antennaindex,
                                            int scansourceindex, int order,
                                            double *delaycoeffs) const
{
  if(antennaindex < 0 || static_cast<size_t>(antennaindex) >= d_tables.size())
    throw std::out_of_range("DelayModel: antenna index out of range");

  const DelayTable &tab = *d_tables[antennaindex];
  int scan = d_basescan + scanindex;
  double t0 = d_startabs + offsettime;

  // DELAY(us) = -delay_sec * 1e6 (SFXC and DiFX use opposite conventions).
  if(order == 0)
  {
    double d = -1e6 * tab.delay(scan, scansourceindex, t0 + timespan / 2.0);
    delaycoeffs[0] = 0.0;
    delaycoeffs[1] = d;
    return true;
  }
  if(order != 1)
    throw std::runtime_error("DelayModel: only interpolator order 0 and 1 "
                             "are supported");

  // Same 3-point linear fit as mpifxcorr Model::calculateDelayInterpolator.
  double d0 = -1e6 * tab.delay(scan, scansourceindex, t0);
  double d1 = -1e6 * tab.delay(scan, scansourceindex, t0 + timespan / 2.0);
  double d2 = -1e6 * tab.delay(scan, scansourceindex, t0 + timespan);

  delaycoeffs[0] = (d2 - d0) / numincrements;
  delaycoeffs[1] = d0 + (d1 - (delaycoeffs[0] * numincrements / 2.0 + d0)) / 3.0;
  return true;
}
