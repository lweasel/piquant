import pandas as pd
import os.path

PREQUANT_RESOURCE_TYPE = "prequant"
QUANT_RESOURCE_TYPE = "quant"

OVERALL_USAGE = "overall_usage"

COMMAND = "command"
REAL_TIME_SECS = "real"
USER_MODE_SECS = "user"
KERNEL_MODE_SECS = "sys"
MAX_RESIDENT_MEMORY = "max_mem"

USAGE_STATS = [REAL_TIME_SECS, USER_MODE_SECS,
               KERNEL_MODE_SECS, MAX_RESIDENT_MEMORY]


def get_usage_summary(usage_file):
    usage_info = pd.read_csv(
        usage_file, header=None,
        names=[COMMAND, REAL_TIME_SECS, USER_MODE_SECS,
               KERNEL_MODE_SECS, MAX_RESIDENT_MEMORY])

    return pd.DataFrame([{
        REAL_TIME_SECS: usage_info[REAL_TIME_SECS].sum(),
        USER_MODE_SECS: usage_info[USER_MODE_SECS].sum(),
        KERNEL_MODE_SECS: usage_info[KERNEL_MODE_SECS].sum(),
        MAX_RESIDENT_MEMORY: usage_info[MAX_RESIDENT_MEMORY].max()
    }])


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
