# automatically generated by the FlatBuffers compiler, do not modify

# namespace: feather_v1_fbs

import flatbuffers
from flatbuffers.compat import import_numpy
np = import_numpy()

class CategoryMetadata(object):
    __slots__ = ['_tab']

    @classmethod
    def GetRootAsCategoryMetadata(cls, buf, offset):
        n = flatbuffers.encode.Get(flatbuffers.packer.uoffset, buf, offset)
        x = CategoryMetadata()
        x.Init(buf, n + offset)
        return x

    # CategoryMetadata
    def Init(self, buf, pos):
        self._tab = flatbuffers.table.Table(buf, pos)

    # The category codes are presumed to be integers that are valid indexes into
    # the levels array
    # CategoryMetadata
    def Levels(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(4))
        if o != 0:
            x = self._tab.Indirect(o + self._tab.Pos)
            from feather_v1_fbs.PrimitiveArray import PrimitiveArray
            obj = PrimitiveArray()
            obj.Init(self._tab.Bytes, x)
            return obj
        return None

    # CategoryMetadata
    def Ordered(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(6))
        if o != 0:
            return bool(self._tab.Get(flatbuffers.number_types.BoolFlags, o + self._tab.Pos))
        return False

def CategoryMetadataStart(builder): builder.StartObject(2)
def CategoryMetadataAddLevels(builder, levels): builder.PrependUOffsetTRelativeSlot(0, flatbuffers.number_types.UOffsetTFlags.py_type(levels), 0)
def CategoryMetadataAddOrdered(builder, ordered): builder.PrependBoolSlot(1, ordered, 0)
def CategoryMetadataEnd(builder): return builder.EndObject()
