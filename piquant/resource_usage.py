import math
import pandas as pd
import os.path

PREQUANT_RESOURCE_TYPE = "prequant_usage"
QUANT_RESOURCE_TYPE = "quant_usage"
OVERALL_USAGE_PREFIX = "overall"


class _ResourceUsageStatistic(object):
    def __init__(self, name, title, format_string, value_extractor):
        self.name = name
        self.title = title
        self.format_string = format_string
        self.value_extractor = value_extractor

    def get_value(self, usage_df):
        return self.value_extractor(usage_df[self.name])

    def stat_range(self, vals_range):
        max_val = math.ceil(vals_range[1] * 2) / 2.0
        return (0, max_val + 0.01)


_RESOURCE_USAGE_STATS = []

_RESOURCE_USAGE_STATS.append(_ResourceUsageStatistic(
    "real-time", "Total elapsed real time", "%e",
    lambda x: math.log10(x.sum())))

_RESOURCE_USAGE_STATS.append(_ResourceUsageStatistic(
    "user-time", "Total user mode time", "%U",
    lambda x: math.log10(x.sum())))

_RESOURCE_USAGE_STATS.append(_ResourceUsageStatistic(
    "sys-time", "Total kernel mode time", "%S",
    lambda x: math.log10(x.sum())))

_RESOURCE_USAGE_STATS.append(_ResourceUsageStatistic(
    "max-memory", "Maximum resident memory", "%M",
    lambda x: x.max() / 1048576.0))


def get_resource_usage_statistics():
    return set(_RESOURCE_USAGE_STATS)


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