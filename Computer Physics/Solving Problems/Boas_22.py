from __future__ import print_function,division
def search(l,number): #l�� ��Ŀ��� �����ִ��� Ȯ��
    for i in range(len(l)):
        if l[i] == number:
            return True
    return False
def count_box(a,N):
    global N_count
    if N == len(Human)-2:
        pass
        print(Human)
    if N == 0:
        for i in range(Box):
            ##
            for j in range(Box):
                if search(Human,j+1) == False:
                    N_count += 1
                    break
            ##
            Human[0] += 1
            
    else:
        for i in range(Box):
            count_box(a,N-1)
            Human[N] += 1
            Human[N-1] = 1
        
N_count = 0
Box = 5
Human = [1,1,1,1,1,1,1,1]
count_box(Human,len(Human)-1)
print("���� Ȯ����",N_count/Box**len(Human))
print("�̷� Ȯ����",((Box-1)**len(Human))/(Box**len(Human)))
    
