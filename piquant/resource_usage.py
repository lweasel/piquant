import math
import pandas as pd
import os.path

PREQUANT_RESOURCE_TYPE = "prequant_usage"
QUANT_RESOURCE_TYPE = "quant_usage"
TIME_USAGE_TYPE = "time"
MEMORY_USAGE_TYPE = "memory"


class _ResourceUsageStatistic(object):
    def __init__(self, name, usage_type, title,
                 format_string, units, value_extractor):
        self.name = name
        self.usage_type = usage_type
        self.title = title
        self.format_string = format_string
        self.units = units
        self.value_extractor = value_extractor

    def get_value(self, usage_df):
        return self.value_extractor(usage_df[self.name])

    def stat_range(self, vals_range):
        max_val = math.ceil(vals_range[1] * 2) / 2.0
        return (0, max_val + 0.01)

    def get_axis_label(self):
        return "{t} ({u})".format(t=self.title, u=self.units)


_RESOURCE_USAGE_STATS = []

_RESOURCE_USAGE_STATS.append(_ResourceUsageStatistic(
    "real-time", TIME_USAGE_TYPE,
    "Log10 total elapsed real time", "%e", "s",
    lambda x: math.log10(x.sum())))

_RESOURCE_USAGE_STATS.append(_ResourceUsageStatistic(
    "user-time", TIME_USAGE_TYPE,
    "Log10 total user mode time", "%U", "s",
    lambda x: math.log10(x.sum())))

_RESOURCE_USAGE_STATS.append(_ResourceUsageStatistic(
    "sys-time", TIME_USAGE_TYPE,
    "Log10 total kernel mode time", "%S", "s",
    lambda x: math.log10(x.sum())))

_RESOURCE_USAGE_STATS.append(_ResourceUsageStatistic(
    "max-memory", MEMORY_USAGE_TYPE,
    "Maximum resident memory", "%M", "Gb",
    lambda x: x.max() / 1048576.0))


def get_resource_usage_statistics():
    return set(_RESOURCE_USAGE_STATS)


def get_time_usage_statistics():
    return [rus for rus in _RESOURCE_USAGE_STATS
            if rus.usage_type == TIME_USAGE_TYPE]


def get_memory_usage_statistics():
    return [rus for rus in _RESOURCE_USAGE_STATS
            if rus.usage_type == MEMORY_USAGE_TYPE]


def get_time_command(resource_type):
    output_file = get_resource_usage_file(resource_type)
    format_string = ",".join(
        [rus.format_string for rus in _RESOURCE_USAGE_STATS])
    return ("/usr/bin/time -f \"\\\"%C\\\",{format_string}\" " +
            "-o {output_file} -a ").format(
        format_string=format_string, output_file=output_file)


def get_usage_summary(usage_file):
    usage_info = pd.read_csv(
        usage_file, header=None,
        names=["command"] + [rus.name for rus in _RESOURCE_USAGE_STATS])

    return pd.DataFrame([
        {rus.name: rus.get_value(usage_info) for rus in _RESOURCE_USAGE_STATS}
    ])


def get_resource_usage_file(resource_type, prefix=None, directory=None):
    file_name = resource_type + ".csv"
    if prefix:
        file_name = prefix + "_" + file_name
    if directory:
        file_name = os.path.join(directory, file_name)
    return file_name


def write_usage_summary(usage_file_name, usage_summary):
    with open(usage_file_name, "w") as out_file:
        usage_summary.to_csv(out_file, index=False)
