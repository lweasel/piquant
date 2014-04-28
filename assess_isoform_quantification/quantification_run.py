import os.path


def _get_reads_name(length, depth, paired_end, errors, bias):
    return "{d}x_{l}b_{p}_{e}_{b}".format(
        d=depth, l=length,
        p="pe" if paired_end else "se",
        e="errors" if errors else "no_errors",
        b="bias" if bias else "no_bias")


def get_run_name(quant_method, length, depth, paired_end, errors, bias):
    return "{qm}_{r}".format(
        qm=quant_method.get_name(),
        r=_get_reads_name(length, depth, paired_end, errors, bias))


def get_reads_dir(output_dir, length, depth, paired_end, error, bias):
    return output_dir + os.path.sep + \
        _get_reads_name(length, depth, paired_end, error, bias)


def get_run_dir(output_dir, quant_method, length, depth,
                paired_end, error, bias):
    return output_dir + os.path.sep + \
        get_run_name(quant_method, length, depth, paired_end, error, bias)
