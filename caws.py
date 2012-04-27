"""Code to process data from the CalTrans Automated Warning System (CAWS).

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import copy
import csv
import datetime
import os
import sys

import matplotlib.pyplot as pyplot
import myplot
import Cdf
import Pmf
import thinkstats

import rpy2.robjects as robjects
r = robjects.r


def run_model(model, print_flag=True):
    """Submits model to r.lm and returns the result."""
    model = r(model)
    res = r.glm(model, family='poisson')
    if print_flag:
        print_summary(res)
    return res


def print_summary(res):
    """Prints results from r.lm (just the parts we want)."""
    flag = False
    lines = r.summary(res)
    lines = str(lines)

    for line in lines.split('\n'):
        # skip everything until we get to coefficients
        if line.startswith('Coefficients'):
            flag = True
        if line.startswith('Signif'):
            continue
        if flag:
            print line
    print


def read_csv(filename, head_lines=1):
    """Reads a CSV file, returns a list of columns.

    filename: string filename
    head_lines: number of header lines to skip
    """
    fp = open(filename)
    reader = csv.reader(fp)

    for i in range(head_lines):
        header = reader.next()

    rows = [t for t in reader]
    fp.close()

    cols = zip(*rows)
    return cols


def write_csv(filename, header, data):
    """Writes a CSV file

    filename: string filename
    header: list of strings
    data: list of rows
    """
    fp = open(filename, 'w')
    writer = csv.writer(fp)
    writer.writerow(header)

    for t in data:
        writer.writerow(t)
    fp.close()


def print_cols(cols):
    """Prints the index and first two elements for each column.

    cols: list of columns
    """
    for i, col in enumerate(cols):
        print i, col[0], col[1]


def investigate_cols(cols):
    data = zip(*cols)
    for row in data:
        if row[7] > 4:
            print row


def filter_cols(cols, indices):
    """Selects columns from a dataset and returns a map from name to column.

    cols: list of columns
    indices: list of (name, index) pairs
    """
    col_dict = {}
    for name, index in indices:
        col_dict[name] = cols[index]
    return col_dict


def make_objects(col_dict, constructor):
    """Turns a column dictionary into a list of objects.

    col_dict: map from attribute name to column
    constructor: function that makes the objects

    Returns: list of objects
    """
    n = len(col_dict.values()[0])

    objs = []
    for i in range(n):
        obj = constructor()
        for name, col in col_dict.iteritems():
            setattr(obj, name, col[i])

        objs.append(obj)

    return objs


class CawsWeather(object):
    def __init__(self):
        self.fog_days = {}

    def read_line(self, reader):
        try:
            return reader.next()
        except StopIteration:
            return None

    def read_csv(self, year):
        filename = 'Weather Data/%d Weather Data.csv' % year
        fp = open(filename)
        reader = csv.reader(fp)
        header = self.read_line(reader)

        count = 0
        month = 0
        while True:
            month += 1
            # read a NOAA sheet and ignore it
            caws_header, lines = self.read_sheet(reader)

            # read a CAWS sheet and parse it
            noaa_header, lines = self.read_sheet(reader)
            count += self.read_month(year, month, caws_header, lines)

            if noaa_header is None:
                break

        fp.close()
        print year, count
        return count

    def read_sheet(self, reader):
        lines = []
        while True:
            t = self.read_line(reader)
            if t is None:
                return t, lines
            
            if len(t) > 2 and t[1] == 'Time':
                return t, lines

            lines.append(t)

    def read_month(self, year, month, header, lines):
        count = 0
        offset = len(header) - 20
        for line in lines:
            try:
                day = int(line[0])
                time = int(line[1])
                weather = line[6+offset]
            except (IndexError, ValueError):
                continue

            count += 1
            key = (year, month, day)
            val = (time, weather)
            self.fog_days.setdefault(key, []).append(val)

        print year, month, count
        return count

    def print_fog_days(self):
        for key, vals in sorted(self.fog_days.iteritems()):
            print key, len(vals)

    def lookup(self, year, month, day):
        if year < 1997 or year > 2001:
            return 'NA'

        key = (year, month, day)
        vals = self.fog_days.get(key, [])
        return vals


class AccidentVehicle(object):
    """Represents a vehicle involved in an accident."""
    
    def dump(self):
        for attr in self.__dict__:
            print attr, getattr(self, attr)


class Accident(object):
    """Encapsulates accident data and provides accessor methods."""

    def __init__(self):
        """Initializes Accident structures.

        subsets: map from subset name to an AV dictionary
        index: map from subset name to list of three Hist objects
               (maps from from (year, month, day) to count)
        """
        all_dict = self.read_accident_data()
        treatment, control = self.partition_av_dict(all_dict)
        self.subsets = dict(treatment=treatment, control=control)

        # index the subsets
        self.index = {}
        for name, subset in self.subsets.iteritems():
            print name
            self.index[name] = self.index_accidents(subset)

    def lookup(self, date, subset_name='control'):
        """Returns the number of accident-vehicles on a given date.

        date: datetime object
        subset_name: string

        Returns a triple of (all accidents, injury accidents, fatal accidents)
        """
        all_hist, injury_hist, fatality_hist = self.index[subset_name]
        key = date.year, date.month, date.day
        return (all_hist.Freq(key), 
                injury_hist.Freq(key), 
                fatality_hist.Freq(key))

    def read_accident_data(self):
        filenames = ['I5before.csv', 'I5after.csv', 
                     'sr120before.csv', 'sr120after.csv']

        all_dict = {}

        for filename in filenames:
            av_dict = self.process_accident_file(filename)
            all_dict.update(av_dict)

        return all_dict

    def process_accident_file(self, filename, data_dir='data'):
        """Reads an accident file.

        Returns a map from date-time string to list of AccidentVehicle
        objects.

        I have adjusted the time of two accidents by 1 minute so that
        the (date, time) keys are unique.
        """
        file_indices = [('rte', 1),
                        ('loc', 3),
                        ('soh', 13),
                        ('date', 15),
                        ('time', 16),
                        ('code', 17),
                        ('num_vehicles', 25),
                        ('pt', 26),
                        ('accidents', 39),
                        ('pdo', 40),
                        ('injury', 41),
                        ('fatality', 42),
                        ]

        filename = os.path.join(data_dir, filename)

        # make a list of columns
        cols = read_csv(filename, head_lines=3)

        # select column and make a map from name to column
        col_dict = filter_cols(cols, file_indices)

        # make the columns into objects
        avs = make_objects(col_dict, AccidentVehicle)

        av_dict = self.make_av_dict(avs)

        total = self.check_av_dict(av_dict)

        return av_dict

    def make_av_dict(self, avs):
        """Groups objects by date time key.

        avs: list of AccidentVehicle objects

        returns: map from (date, time) string to list of AccidentVehicle objects
        """
        av_dict = {}
        for av in avs:
            if av.time:
                key = '%s %4.4d' % (av.date, int(av.time))
                av_dict.setdefault(key, []).append(av)
        return av_dict

    def check_av_dict(self, av_dict):
        """Counts the total number of accidents in an av_dict.

        av_dict: map from (date, time) string to list of AccidentVehicles
        """
        total = 0

        for key, avs in av_dict.iteritems():
            #print key, len(avs)
            for av in avs:
                total += int_or_empty(av.accidents)
                #print av.loc, av.pt, av.num_accidents

        return total

    def partition_av_dict(self, av_dict):
        """Partitions an av_dict into control and treatment areas.

        Control area is 5N and 120E.  Treatment area is 5S and 120W.

        Returns two av_dicts.
        """
        def get_rte(avs):
            av = avs[0]
            rte = av.rte + av.soh
            return rte

        treatment = {}
        control = {}

        for key, avs in av_dict.iteritems():
            rte = get_rte(avs)
            if rte in ['5S', '120W']:
                treatment[key] = avs
            elif rte in ['5N', '120E']:
                control[key] = avs
            else:
                raise KeyError('unknown highway')

        return treatment, control

    def plot_data(self, root='caws.accident'):
        """Plots a time series of monthly accidents.

        root: string prefix of the output files.
        """
        pyplot.clf()
        for name, av_dict in self.subsets.iteritems():
            hist = self.count_accidents(av_dict)
            
            print name, 'Total accidents', hist.Total()
            years, counts = zip(*sorted(hist.Items()))
            pyplot.plot(years, counts, label=name)

        myplot.Save(root=root,
                    title='Monthly Accident Counts',
                    xlabel='Year',
                    ylabel='Number of accidents',
                    axis=[1991.5, 2002.5, 0, 40])

    def count_accidents(self, av_dict):
        """Counts the number of accidents in each month.

        Returns a Hist object that maps months (as floating-point
        years) to number of accidents.
        """
        hist = Pmf.Hist()
        for key, avs in av_dict.iteritems():
            year, month, day = convert_date_key(key)
            year_frac = convert_year_month(year, month)
            hist.Incr(year_frac)
        return hist

    def plot_agg_data(self, traffic, index=0, root='caws.accident'):
        """Plots a time series of monthly accidents.

        index: which histogram to get data from
            0 = all accidents
            1 = injury or fatal
            2 = fatal
        root: string prefix of the output files.
        """
        # upper bounds for the different plots
        d = {0:10, 1:3, 2:0.25}
        ymax = d[index]

        pyplot.clf()

        # draw lines for significant dates
        # 3-25-1996: speed limit increased on I-5
        # 4-22-1996: speed limit increased on SR-120
        # 11-?-1996: caws turned on
        for date in [(1996, 3.8), (1996, 4.75), (1996, 11.5)]:
            year, month = date
            year_frac = convert_year_month(year, month)
            xs = [year_frac, year_frac]
            ys = [0, ymax]
            pyplot.plot(xs, ys, color='red', linewidth=2, alpha=0.2)

        style = [dict(color='blue',
                      markerfacecolor=(0, 0, 1, 0.5), 
                      markeredgecolor=(0, 0, 1, 0.5)),
                 dict(color='green', 
                      markerfacecolor=(0, 0.5, 0, 0.5),
                      markeredgecolor=(0, 0.5, 0, 0.5),
                      )
                 ]

        # draw data lines
        i = 0
        for name, hists in self.index.iteritems():
            hist = self.aggregate_accidents(traffic, hists[index])            
            years, counts = zip(*sorted(hist.Items()))
            pyplot.plot(years, counts, 
                        label=name, 
                        linewidth=2, 
                        alpha=0.5,
                        marker='o',
                        markersize=3,
                        **style[i])
            i += 1

        myplot.Save(root=root,
                    title='Quarterly Accident Rates',
                    xlabel='Year',
                    ylabel='Accidents per million cars',
                    axis=[1991.5, 2002.5, 0, ymax])

    def aggregate_accidents(self, traffic, hist):
        """Counts the number of accidents in each month.

        Returns a Hist object that maps months (as floating-point
        years) to number of accidents.
        """
        agg = self.empty_hist(1992.25, 2002.0)
        for key, count in hist.Items():
            year, month, day = key
            
            # estimate the traffic for this quarter
            aadt = traffic.lookup(year) * 365 / 4.0

            # round off the month part to the nearest quarter
            # each data point reflects the PREVIOUS three months
            month = lump_month(month, 3)
            year_frac = convert_year_month(year, month) + 0.25

            # convert count to a rate per million cars
            rate = float(count) / aadt * 1000000
            agg.Incr(year_frac, rate)
        return agg

    def empty_hist(self, year1, year2, step=0.25):
        """Makes a Hist with zero entries for each quarter.

        year1: first quarter (float year)
        year2: last quarter, inclusive
        step: step size, float years

        Returns Hist object
        """
        hist = Pmf.Hist()
        year = year1
        while year <= year2:
            hist.Set(year, 0)
            year += step
        return hist

    def index_accidents(self, av_dict):
        """Makes a map from (year, month, day) to number of accidents.

        Returns three Hist objects: all accidents, injury or fatal, fatal. 
        Each one maps (year, month, day) to number of accidents.
        """
        all_hist = Pmf.Hist()
        injury_hist = Pmf.Hist()
        fatality_hist = Pmf.Hist()

        for key, avs in av_dict.iteritems():
            date = convert_date_key(key)

            # all_hist includes injury and fatal.
            # injury_hist accidents includes fatal.
            all_hist.Incr(date)
            injury = av_attr(avs, 'injury')
            fatality = av_attr(avs, 'fatality')
            if injury or fatality:
                injury_hist.Incr(date)
            if fatality:
                fatality_hist.Incr(date)

        print all_hist.Total(), injury_hist.Total(), fatality_hist.Total()
        return all_hist, injury_hist, fatality_hist


def lump_month(month, f):
    return 1 + int((month-1)/f) * f

def av_attr(avs, attr):
    """Gets an attribute from the first AccidentVehicle in a list.

    avs: list of AccidentVehicle object.
    attr: string attr

    Returns an integer.  The empty string is considered 0.
    """
    return int_or_empty(getattr(avs[0], attr))


def convert_date_key(key):
    """Converts a date-time string into a (year, month, day) triple of int.

    key: string containing a date and time, separated by a space
    """
    date, time = key.split()
    month, day, year = date.split('/')
    return int(year), int(month), int(day)


def convert_year_month(year, month):
    """Converts year and month to a floating-point year with a fraction."""
    year_frac = float(year) + (month-1)/12.0
    return year_frac


def int_or_empty(s):
    """Converts a string to integer.

    The empty string is considered 0.
    """
    if s == '':
        return 0
    return int(s)


class TrafficLoc(object):
    """Represents a record of average traffic at a location."""


class Traffic(object):
    """Encapsulates traffic data from CalTrans.

    Data files downloaded from http://traffic-counts.dot.ca.gov/
    on April 12, 2012.  Converted xls to csv.

    Traffic reports are AADT: annual average daily traffic; that
    is, an estimate of daily traffic in both directions, averaged
    over the year.

    For each location there are four reported numbers:

    ahead: measurement from the sensor positioned before the
           nominal measurement point (normally to the N or E)

    back: measurement from the sensor positioned after the
           nominal measurement point (normally to the S or W)

    ahead_peak: estimated daily average for peak month, ahead position
    back_peak: estimated daily average for peak month, back position
    """
    # range of years
    years = range(1992, 2003)
    
    # map from mile marker to label
    locs = {'3.32':'SR-120', '14.83':'merged', '17.52':'I-5'}

    def __init__(self):
        self.all_dict = self.read_traffic_data()

    def lookup(self, year, loc='14.83', attr='back'):
        """Looks up the AADT for a year and location.

        year: int year
        loc: string, one of the keys in locs
        attr: which variable to look up

        Returns int cars per day.
        """
        for aadt in self.all_dict[year]:
            if aadt.loc == loc:
                return int(getattr(aadt, attr))

    def read_traffic_data(self):
        """Read the traffic data files.

        Returns a map from year to list of TrafficLoc objects.
        """
        all_dict = {}
        for year in self.years:
            aadt_dict = self.process_aadt_file(year)
            all_dict.update(aadt_dict)
        return all_dict

    def process_aadt_file(self, year, data_dir = 'data'):
        """Reads a csv file with AADT data.

        year: int year
        data_dir: string name of directory

        Returns a map from year to list of TrafficLoc objects.        
        """
        filename = '%dadt.csv' % year
        filename = os.path.join(data_dir, filename)

        file_indices = [('rte', 1),
                        ('loc', 4),
                        ('desc', 5),
                        ('back_peak', 7),
                        ('back', 8),
                        ('ahead_peak', 10),
                        ('ahead', 11),
                        ]

        # make a list of columns
        cols = read_csv(filename)

        # select column and make a map from name to column
        col_dict = filter_cols(cols, file_indices)

        # make the columns into objects
        aadts = make_objects(col_dict, TrafficLoc)

        aadt_dict = {}
        for aadt in aadts:
            if aadt.rte in ['5', '120']:
                if aadt.loc in self.locs:
                    aadt_dict.setdefault(year, []).append(aadt)

        return aadt_dict

    def plot_data(self, root='caws.traffic'):
        """Makes a plot of AADT for each location."""
        pyplot.clf()
        series = {}

        for loc, name in self.locs.iteritems():
            for year in self.years:
                adt = self.lookup(year, loc) / 1000
                series.setdefault(name, []).append(adt)

        # TODO: fix the year labels
        for name, adts in series.iteritems():
            pyplot.plot(self.years, adts, label=name)

        myplot.Save(root=root,
                    title='Traffic volume',
                    xlabel='Year',
                    ylabel='AADT',
                    axis=[1991.5, 2002.5, 0, 160])


class WeatherDay(object):
    """Represents weather data for one day."""


class Weather(object):
    """Encapsulates weather data."""

    def __init__(self):
        self.wds = self.read_weather_data()

    def read_weather_data(self, 
                          filename='ksck_weather_29835.csv',
                          data_dir = 'data'):
        """Reads the weather file downloaded from NOAA:

        http://www.ncdc.noaa.gov/cdo-web/search;
        jsessionid=6EFC504444861FF6F540A5B16456A0D1.lwf1#t=secondTabLink

        KSCK: Stockton airport.
        Daily

        Start Date: 1/1/1993
        End Date:	 12/31/2002
 	Requested Data
Stations:	GHCND:USW00023237 - STOCKTON METROPOLITAN AIRPORT
Data Types:	WT18 - Snow, snow pellets, snow grains, or ice crystals
WT05 - Hail (may include small hail)
WT19 - Unknown source of precipitation 
WT03 - Thunder
AWND - Average daily wind speed (tenths of meters per second)
WT04 - Ice pellets, sleet, snow pellets, or small hail" 
WT09 - Blowing or drifting snow
WT14 - Drizzle
WT16 - Rain (may include freezing rain, drizzle, and freezing drizzle)" 
WT07 - Dust, volcanic ash, blowing dust, blowing sand, or blowing obstruction
WT08 - Smoke or haze 
WT22 - Ice fog or freezing fog
TSUN - Daily total sunshine (minutes)
PRCP - Precipitation (tenths of mm)
WT21 - Ground fog 
WT13 - Mist
WT02 - Heavy fog or heaving freezing fog (not always distinguished from fog)
WT01 - Fog, ice fog, or freezing fog (may include heavy fog)
WT10 - Tornado, waterspout, or funnel cloud" 
        """
        file_indices = [('date', 1),
                        ('prcp', 2),     # tenths of mm
                        ('sun', 3),      # sunshine in minutes
                        ('wind', 4),     # AWND
                        ('fog', 5),      # WT01
                        ('hfog', 6),     # WT02
                        ]

        filename = os.path.join(data_dir, filename)

        # make a list of columns
        cols = read_csv(filename, head_lines=1)

        # select column and make a map from name to column
        col_dict = filter_cols(cols, file_indices)

        # make the columns into objects
        wds = make_objects(col_dict, WeatherDay)
        for wd in wds:
            self.clean_wd(wd)
        return wds

    def clean_wd(self, wd):
        """Converts the fog attributes to int.

        Rewrites attributes fog and hfog

        wd: WeatherDay object
        """
        d = {'1':1, '9999':0}
        wd.fog = d[wd.fog]
        wd.hfog = d[wd.hfog]
        wd.prcp = float(wd.prcp) / 10.0

    def get_wds(self):
        """Generator that iterates WeatherDays.

        Yields (datetime, WeatherDay) pairs
        """
        format = '%Y%m%d'
        for wd in self.wds:
            date = datetime.datetime.strptime(wd.date, format)
            yield date, wd

    def count_days(self, attr, thresh):
        """Counts the number of days of fog/precipitation in each month.

        attr: string attribute to count
        thresh: value to exceed

        Returns: Hist object that maps years to count.
        """
        hist = Pmf.Hist()
        format = '%Y%m%d'
        for wd in self.wds:
            date = datetime.datetime.strptime(wd.date, format)
            key = date.year, date.month
            if getattr(wd, attr) > thresh:
                hist.Incr(key)

        return hist

    def plot_data(self, root='caws.weather'):
        """Plots number of fog days per month."""
        pyplot.clf()

        def time_series(hist):
            """Convert Hist object to a time series.
            
            Returns a list of (year, count) pairs.
            """
            data = []
            for key, count in sorted(hist.Items()):
                year, month = key
                year_frac = convert_year_month(year, month)
                data.append((year_frac, count))
            return data

        attrs = [('fog', 0), ('hfog', 0), ('prcp', 3)]

        for attr, thresh in attrs:
            data = self.count_days(attr, thresh)
            years, ys = zip(*time_series(data))
            pyplot.plot(years, ys, label=attr)

        myplot.Save(root=root,
                    title='Monthly Days of Weather',
                    xlabel='Year',
                    ylabel='Number of days',
                    axis=[1991.5, 2002.5, 0, 35])


class Merger(object):
    """Contains methods for writing and reading the merged data sets."""

    header = ['year', 
              'month', 
              'day', 
              'i', 
              'fog',
              'hfog',
              'precip', 
              'adt',
              'accidents',
              'injuries',
              'fatalities',
              'accidents_per',
              'injuries_per',
              'fatalities_per',
              'speed',
              'caws'
              ]

    def create_merged_data(self, plot_flag=False):
        """Reads raw data files and writes a merged dataset for analysis."""
        accident = Accident()
        traffic = Traffic()
        weather = Weather()

        if plot_flag:
            make_plots(accident, traffic, weather)

        header, data = self.merge_data(accident, traffic, weather, 'control')
        filename='data/control_data.csv'
        print 'Writing', filename
        write_csv(filename, header, data)

        header, data = self.merge_data(accident, traffic, weather, 'treatment')
        filename='data/treatment_data.csv'
        print 'Writing', filename
        write_csv(filename, header, data)


    def merge_data(self, accident, traffic, weather, subset_name,
                   loc='14.83', end=(2002,3,31)):
        """Generates a table of data for analyis.

        accident: Accident object with accident data
        traffic: Traffic object
        weather: Weather object
        subset_name: string
        loc: mile marked to use for traffic volume
        end: (year, month, day) triple for the stop date


        loc = '14.83' selects total traffic after the merge.
        end is the last day of the last month where we have accident data.
        """
        data = []
        for i, (date, wd) in enumerate(weather.get_wds()):
            year, month, day = date.year, date.month, date.day

            adt = traffic.lookup(year, loc)

            accidents, injuries, fatalities = accident.lookup(date, subset_name)
            accidents_per = accidents * 1000000.0 / adt
            injuries_per = injuries * 1000000.0 / adt
            fatalities_per = fatalities * 1000000.0 / adt

            # more than 3mm of rain
            precip = 1 if int(wd.prcp) >= 3 else 0

            speed_date = datetime.datetime(1996, 1, 15)
            caws_date = datetime.datetime(1996, 11, 15)

            speed = 1 if date >= speed_date else 0
            caws = 1 if date >= caws_date else 0

            t = [year, 
                 month, 
                 day,
                 i, 
                 wd.fog,
                 wd.hfog,
                 wd.prcp,
                 adt/1000.0, 
                 accidents, 
                 injuries, 
                 fatalities, 
                 accidents_per, 
                 injuries_per, 
                 fatalities_per, 
                 speed, 
                 caws
                 ]
            data.append(t)

            if (year, month, day) == end:
                break

        return self.header, data


    def read_merged_file(self, filename, data_dir='data'):
        """Reads the merged dataset.

        Returns a map from attribute name to column of data.
        """
        file_indices = [(name, i) for i, name in enumerate(self.header)]
        filename = os.path.join(data_dir, filename)

        # make a list of columns
        cols = read_csv(filename, head_lines=1)
        cols = [float_col(col) for col in cols]

        #investigate_cols(cols)

        # select column and make a map from name to column
        col_dict = filter_cols(cols, file_indices)

        return col_dict


def int_col(col):
    """Converts a column of data to int."""
    return [int(x) for x in col]


def float_col(col):
    """Converts a column of data to float."""
    return [float(x) for x in col]


def inject_col_dict(col_dict, prefix=''):
    """Copies a data columns into the R global environment.
    
    col_dict: map from attribute name to column of data
    prefix: string prepended to the attribute names
    """
    for name, col in col_dict.iteritems():
        robjects.globalEnv[prefix+name] = robjects.FloatVector(col)


def split_col_dict(col_dict, index):
    """Splits a column dictionary into a before and after.
    
    col_dict: map from attribute name to column of data
    index: int index where the split should go

    Returns (before, after) a pair of col_dicts.
    """
    before = copy.deepcopy(col_dict)
    after = copy.deepcopy(col_dict)

    for key, col in before.iteritems():
        before[key] = col[:index]
    
    for key, col in after.iteritems():
        after[key] = col[index:]
    
    return before, after


def accident_cdf(col_dict, attr='accidents'):
    """Makes a CDF of the number of accidents per day.

    col_dict: map from attribute name to column of data
    attr: string attribute name, one of accidents, injury, fatal
    """
    accidents = col_dict[attr]
    cdf = Cdf.MakeCdfFromList(accidents)
    print thinkstats.MeanVar(accidents)
    return cdf


def plot_accident_cdfs():
    """Plots CDF of accident counts for the control and treatment areas.

    Before and after the date CAWS was deployed (known to be November 1996).
    """
    cdfs = []
    for label in ['control', 'treatment']:
        print label
        filename = label + '_data.csv'
        col_dict = process_merged_file(filename)

        # November 15, 1996
        before, after = split_col_dict(col_dict, 1780)
        print 'before'
        cdf = accident_cdf(before, 'accidents')
        cdf.name = label + ' before'
        cdfs.append(cdf)

        print 'after'
        cdf = accident_cdf(after, 'accidents')
        cdf.name = label + ' after'
        cdfs.append(cdf)

    myplot.Cdfs(cdfs, 
                root='caws.poisson',
                transform='exponential',
                title='CCDF of Accident Counts',
                xlabel='Number of accidents',
                ylabel='Complementary CDF')


def make_plots(accident, traffic, weather):
    """Plots time series data."""
    weather.plot_data()
    accident.plot_data()
    for index in range(3):
        root = 'caws.accident.%d' % index
        accident.plot_agg_data(traffic, index, root)
    traffic.plot_data()


def run_model_on_partitions(model):
    """Runs models on the control and treatment subsets, before and after.

    model: string representation of the model, in R format.

    Prints the results
    """
    merger = Merger()

    for label in ['control', 'treatment']:
        print label, model
        filename = label + '_data.csv'
        col_dict = merger.read_merged_file(filename)
        inject_and_run_model(model, col_dict)
        before, after = split_col_dict(col_dict, 1414)
        print label, 'before:', model
        inject_and_run_model(model, before)
        print label, 'after:', model
        inject_and_run_model(model, after)
        print ''


def inject_and_run_model(model, col_dict):
    """Runs regression models in R."""
    inject_col_dict(col_dict)
    run_model(model)


def main(script):
    cw = CawsWeather()
    # TODO: read the 2002 data, which has just CAWS data, only one month
    # of NOAA data
    for year in range(1997, 2002):
        cw.read_csv(year)
    cw.print_fog_days()
    return

    #plot_accident_cdfs()
    #merger = Merger()
    #merger.create_merged_data(plot_flag=True)

    for dep in ['accidents', 'injuries', 'fatalities']:
        model = '%s ~ adt + precip' % dep
        run_model_on_partitions(model)

if __name__ == '__main__':
    main(*sys.argv)

