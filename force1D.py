# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from feon.sa import *

if __name__ == "__main__":
    E = 210e6
    A = 0.005
    I = 5e-5

    n0 = Node(0, 0)
    n1 = Node(2, 0)
    n2 = Node(4, 0)
    n3 = Node(6, 0)
    n4 = Node(8, 0)
    n5 = Node(10, 0)

    e0 = Beam1D11((n0, n1), E, A, I)
    e1 = Beam1D11((n1, n2), E, A, I)
    e2 = Beam1D11((n2, n3), E, A, I)
    e3 = Beam1D11((n3, n4), E, A, I)
    e4 = Beam1D11((n4, n5), E, A, I)


    s = System()
    s.add_nodes(n0, n1, n2, n3, n4,n5)
    s.add_elements(e0, e1, e2, e3,e4)


    s.add_fixed_sup(0) #n0为固定端
    s.add_element_load(0, "Q", -7)
    s.add_element_load(1, "Q", -7)
    s.add_element_load(2, "Q", -7)
    s.add_element_load(3, "Q", -7)
    s.add_element_load(4, "Q",-10)
    s.solve()

    from feon.sa.draw2d import *

    for el in [e0, e1, e2, e3,e4]:
        draw_bar_info(el)






