import numpy as np
import tkinter as tk
from tkinter import filedialog
from openpyxl import load_workbook
def click():
    fname = filedialog.askopenfilename(filetypes=(("Execel files", "*.xlsx"),("all files","*.*")))
    if fname:
        run(fname)
        window.destroy()
def run(dataname):
    load_wb = load_workbook(dataname,data_only=True)
    load_ws = load_wb['References']

    all_values = []
    for row in load_ws.rows: # 엑셀 모든 데이터 불러오기
        row_value = []
        for cell in row:
            row_value.append(cell.value)
        all_values.append(row_value)
    Data = np.array(all_values) #계산하기 위해 Numpy로 데이터 치환
    set_point = 0
    for data_l in Data[0]:
        if data_l == None: #데이터와 데이터 자료 기준점 찾기
            break
        set_point += 1
    for i in range(np.shape(Data)[0]-1):
        try:
            np.savetxt('Data\\'+str(i+1)+'_'+Data[i+1,set_point+2]+'.csv',np.vstack((np.flip(Data[0,1:set_point]),np.flip(Data[i+1,1:set_point]))),delimiter=',')
            print(i,"번째 자료 완성"," 진행도:",100*i/(np.shape(Data)[0]-1),"%")
        except:
            print(i,"  No Data!!")
window=tk.Tk()
window.geometry("200x100")
window.resizable(False,False)
window.title("Converter by KKS")
btn = tk.Button(window,text = "불러오기",command = click)
btn.grid(row=0,column=0)
window=tk.mainloop()