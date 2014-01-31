from schema import Schema, Use


def validate_file_option(file_option, msg):
    msg = "{msg} '{file}'.".format(msg=msg, file=file_option)
    return Schema(open, error=msg).validate(file_option)


def validate_dict_option(dict_option, values_dict, msg):
    msg = "{msg} '{opt}'.".format(msg=msg, opt=dict_option)
    return Schema(Use(lambda x: values_dict[x]), error=msg).\
        validate(dict_option)
