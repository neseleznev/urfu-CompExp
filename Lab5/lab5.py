#Thanks to Roman Dubinin http://github.com/RomanDubinin
#


from math import e
import matplotlib.pyplot as plt
import numpy as np
import sys
from tkinter import Tk, W, Frame, Button, Entry, StringVar, HORIZONTAL, Label, LEFT, CENTER, FLAT, GROOVE, filedialog, IntVar, Checkbutton

from lab5meths import explicit_euler, recount_euler, runge_kutta, get_next_M

methods = [
    (explicit_euler, 'r8-'),
    (recount_euler, 'rs-'),
    (runge_kutta, 'yd-')
]


N = 14
l = 2 + 0.1 * N

def F(x, y, d_y):
    return y + l*x*(1-x) + 2*l + 2

#y(0) = 0; y'(1) = e-1/e+2.1
def real_func(x):
    return 2.1*x*x - 2.1*x + 1*e**x + e**(-x) - 2

def real_func_method(x_points):
    y_points = []
    for i in range(len(x_points)):
        y_points.append(real_func(x_points[i]))

    return y_points

def draw_all(n, meth):
    print(meth)
    y_0 = 0
    h = 1/n
    right_boundary = e - 1/e + l
    plt.close()
    x_points = np.linspace(0, 1, n+1, endpoint=True)
    for i, pair in enumerate(methods):
	    print(i, pair)
	    if meth[i] == 1:
	        method, color = pair
	        last_M = 0
	        y_points = [method(F, x_points, y_0, last_M, h)]
	        next_M = get_next_M(x_points, F, method, y_0, right_boundary, last_M, h)
	        while abs(next_M - last_M) > 0.01:
	            y_points=method(F, x_points, y_0, next_M, h)
	            #print(y_points)
	            last_M = next_M
	            next_M = get_next_M(x_points, F, method, y_0, right_boundary, last_M, h)
	        plt.plot(x_points, y_points, color, label=method.__doc__)
    y_points = real_func_method(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, 'b*-', label='Original')

    legend = plt.legend(loc='upper left', shadow=True, fontsize='x-large')
    legend.get_frame().set_facecolor('#00FFCC')

    mng = plt.get_current_fig_manager()

    mng.window.state('zoomed')
    # mng.full_screen_toggle()
    plt.show()



window = Tk()
# window.iconbitmap('Nik.ico')
window.title('lab')
window.resizable(width=False, height=False)
color1 = 'white'

frame = Frame(window, bg=color1)
frame.pack()

#l1 = Label(frame, text="Лабораторная работа №1.\n1 вариант.",
 #          justify=CENTER, font=("Helvetica", 12), bd=10, bg=color1).grid(row=0, sticky=W)

eul_exp = IntVar()
Checkbutton(frame, text="Эйлер явный", variable=eul_exp).grid(row=1, sticky=W)
eul_rec = IntVar()
Checkbutton(frame, text="Эйлер с пересчетом", variable=eul_rec).grid(row=2, sticky=W)
r_k = IntVar()
Checkbutton(frame, text="Рунге-Кутта", variable=r_k).grid(row=5, sticky=W)


#l2 = Label(frame, text="\nВыберите количество точек:", justify=CENTER,
                 #font=("Arial", 12), bd=0, bg=color1).grid(row=9, sticky=W)
Label(frame, text="\nВыберите количество точек:", justify=CENTER).grid(row=9, sticky=W)

e1 = Entry(frame)
e1.grid(row=9, column=1, sticky=W)
print(e1.get())

button = Button(frame, text="Нарисовать график", font=("Helvetica", 12),
         bg=color1, command=lambda: draw_all(int(e1.get()), [eul_exp.get(), eul_rec.get(), r_k.get()]
)).grid(row=11, sticky=W)

window.mainloop()
