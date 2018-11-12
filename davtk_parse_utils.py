import argparse

# subclass ArgumentParser to throw errors instead of exiting
class ArgumentParserError(Exception):
    pass

class ThrowingArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ArgumentParserError(message)
    def exit(self):
        raise ArgumentParserError("help")
