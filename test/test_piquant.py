import os.path
import piquant.piquant as piq
import piquant.piquant_options as po


def test_get_parameters_dir_returns_correct_path():
    output_dir = "dummy"
    options = {
        po.OUTPUT_DIRECTORY: output_dir
    }

    assert piq._get_parameters_dir(
        options, read_depth=30, read_length=50,
        paired_end=True, bias=False) == \
        output_dir + os.path.sep + "30x_50b_pe_no_bias"
