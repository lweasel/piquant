from schema import Schema


def validate_file_option(file_option, msg):
    msg = "{msg} '{file}'.".format(msg=msg, file=file_option)
    return Schema(open, error=msg).validate(file_option)
