# automatically generated by the FlatBuffers compiler, do not modify

# namespace: feather_v1_fbs

class Encoding(object):
    PLAIN = 0
    # Data is stored dictionary-encoded
    # dictionary size: <INT32 Dictionary size>
    # dictionary data: <TYPE primitive array>
    # dictionary index: <INT32 primitive array>
    #
    # TODO: do we care about storing the index values in a smaller typeclass
    DICTIONARY = 1

