def spaces_to_underscore(string):
    return string.replace(' ', '_')


def get_order_string(ascending):
    return "asc" if ascending else "desc"
