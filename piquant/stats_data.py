import itertools

from . import piquant_options as po


class OptionValuesSets(object):
    def __init__(self, stats_df):
        self.values = {o: stats_df[o.name].value_counts().index.tolist()
                       for o in po.get_multiple_quant_run_options()}

    def _is_degenerate_option(self, option):
        return len(self.values[option]) <= 1

    def _remove_from(self, options, to_remove):
        get_pset = lambda x: x if isinstance(x, set) \
            else (set(x) if isinstance(x, list) else set([x]))
        return get_pset(options) - get_pset(to_remove)

    def _get_fixed_options(self, non_fixed_options):
        fixed_options = self._remove_from(
            po.get_multiple_quant_run_options(), non_fixed_options)
        non_deg_fixed_options = \
            [o for o in fixed_options if not self._is_degenerate_option(o)]
        value_sets = [v for v in itertools.product(
                      *[self.values[o] for o in non_deg_fixed_options])]
        return non_deg_fixed_options, value_sets

    def get_non_degenerate_options(
            self, numeric_only=False, opts_to_remove=None):

        options = po.get_numerical_mqr_options() if numeric_only \
            else po.get_multiple_quant_run_options()
        if opts_to_remove:
            options = self._remove_from(options, opts_to_remove)
        return [o for o in options if not self._is_degenerate_option(o)]

    def exec_for_fixed_option_values_sets(self, func, non_fixed_options, data):
        fixed_options, fo_values_sets = \
            self._get_fixed_options(non_fixed_options)

        for fo_values_set in fo_values_sets:
            fixed_option_values = {}
            data_subset = data
            for i, fixed_o in enumerate(fixed_options):
                fo_value = fo_values_set[i]
                data_subset = data_subset[data_subset[fixed_o.name] == fo_value]
                fixed_option_values[fixed_o] = fo_value

            func(data_subset, fixed_option_values)
