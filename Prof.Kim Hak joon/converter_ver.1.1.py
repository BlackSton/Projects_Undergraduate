import numpy as np
import tkinter as tk
from tkinter import filedialog
def click():
    fname = filedialog.askopenfilenames(filetypes=(("CSV files", "*.csv"),("all files","*.*")))
    for dir in fname:
        run(dir)
def run(dataname):
    if dataname[-3:] == 'txt': # 파일의 확장자 파악
        deli = ' '
    elif dataname[-3:] == 'csv':
        deli = ','
    filename = dataname[location_file(dataname):-4]
    print("open ",filename)
    data = np.genfromtxt(dataname,skip_header=19,skip_footer =34,dtype = 'f',delimiter=deli)
    t = np.zeros((np.shape(data)[1],np.shape(data)[0]),dtype = 'f')
    i = 0
    for line in data:
        j = 0
        for g in line:
            t[j,i] = g
            j = j +1
        i = i + 1
    for i in range(np.shape(data)[1]-1):
        get_data = np.round(np.vstack((t[0],t[i+1])),3)
        print(get_data)
        np.savetxt('Data\\'+filename+'_'+str(i+1)+'.csv',get_data,fmt = '%.8g',delimiter=',')

def location_file(filedata):
    for i in range(len(filedata)):
        if filedata[-(i+1)] == '/':
            return -(i+1)
            
window=tk.Tk()
window.geometry("280x50")
window.resizable(False,False)
window.title("Converter by KKS")
btn = tk.Button(window,text = "불러오기",command = click)
btn.grid(row=0,column=0)
window=tk.mainloop()
np.genfromtxt()