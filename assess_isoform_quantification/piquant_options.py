import quantifiers


def check_quantification_method(data):
    available_methods = quantifiers.get_quantification_methods()
    return available_methods[data]
