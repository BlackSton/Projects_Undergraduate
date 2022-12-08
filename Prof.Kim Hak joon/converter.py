import tkinter as tk
from tkinter import filedialog
import numpy as np
def click():
    fname = filedialog.askopenfilename(filetypes=(("CSV files", "*.csv"),("TXT files", "*.txt")))
    return fname
    if fname:
        run(fname)
        window.destroy()
def run(dataname):
    print(dataname,tk.Entry.get(txt))
    data = np.loadtxt(dataname,delimiter=tk.Entry.get(txt))
    t = np.zeros((np.shape(data)[1],np.shape(data)[0]))
    i = 0
    for line in data:
        j = 0
        for g in line:
            t[j,i] = g
            j = j +1
        i = i + 1
    for i in range(np.shape(data)[1]-1):
        np.savetxt('Data\\Data'+str(i+1)+'.csv',np.vstack((t[0],t[i+1])),delimiter=',')
window=tk.Tk()
window.geometry("280x50")
window.resizable(False,False)
window.title("Converter by KKS")
lbl = tk.Label(window,text = "Delimiter")
lbl.grid(row=0,column=0)
txt = tk.Entry(window)
txt.grid(row=0,column=1)
btn = tk.Button(window,text = "불러오기",command = click)
btn.grid(row=0,column=2)
window=tk.mainloop()