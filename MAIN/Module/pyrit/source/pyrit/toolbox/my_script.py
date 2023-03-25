# coding=utf-8
"""
Testing stuff with networkx
"""

class A:

    def __init__(self):
        self._a = None

    @property
    def a(self):
        return self._a


class B(A):

    def __init__(self):
        super().__init__()

    # @a.setter
    # def a(self, a):
    #     self._a = a


objA = A()

objB = B()

# objB.a = 2

print(objB.a)