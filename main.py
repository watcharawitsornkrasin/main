#######################################################################
#    Python MRI Quality Controlling
#    Created by : Watcharawit Sornkrasin
#    Last Updated : 22 Apr 2021
#######################################################################
from pydicom import dcmread
import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()
filename_IM0001 = open("DICOM\IM_0001","rb")
filename_IM0005 = open("DICOM\IM_0005","rb")
filename_IM0007 = open("DICOM\IM_0007","rb")
filename_IM0008 = open("DICOM\IM_0008","rb")
filename_IM0009 = open("DICOM\IM_0009","rb")
filename_IM0010 = open("DICOM\IM_0010","rb")
filename_IM0011 = open("DICOM\IM_0011","rb")
fileread_IM0001 = dcmread(filename_IM0001)
fileread_IM0005 = dcmread(filename_IM0005)
fileread_IM0007 = dcmread(filename_IM0007)
fileread_IM0008 = dcmread(filename_IM0008)
fileread_IM0009 = dcmread(filename_IM0009)
fileread_IM0010 = dcmread(filename_IM0010)
fileread_IM0011 = dcmread(filename_IM0011)

pix_array_IM_0001 = fileread_IM0001.pixel_array
pix_array_IM_0005 = fileread_IM0005.pixel_array
pix_array_IM_0007 = fileread_IM0007.pixel_array
pix_array_IM_0008 = fileread_IM0008.pixel_array
pix_array_IM_0009 = fileread_IM0009.pixel_array
pix_array_IM_0010 = fileread_IM0010.pixel_array
pix_array_IM_0011 = fileread_IM0011.pixel_array

#### Processing Geometry Accuracy Test
print('*********  MRI QUALITY CONTROLLING : Geometry Accuracy Test  *********')
print("")
print("Testing from DICOM\IM_0005")
Scale_x_test1 = fileread_IM0005.PixelSpacing[0]
Scale_y_test1 = fileread_IM0005.PixelSpacing[1]
Real_x_test1 = int(fileread_IM0005.Rows)*float(Scale_x_test1)
Real_y_test1 = int(fileread_IM0005.Columns)*float(Scale_y_test1)
print('IM_0005 Resolution = {0:} x {1:}'.format(fileread_IM0005.Rows,fileread_IM0005.Columns))
print('IM_0005 Image size = ',Real_x_test1,'mm','x',Real_y_test1,'mm')
WindowWidth_test1 = float(fileread_IM0005.WindowWidth)
WindowCenter_test1 = float(fileread_IM0005.WindowCenter)
print('IM_0005 Window Width = ',WindowWidth_test1)
print('IM_0005 Window Center = ',WindowCenter_test1)
print()

k_i = 0
k_j = 0
WindowCenter = np.mean(pix_array_IM_0005)
pix_array_test1 = np.copy(pix_array_IM_0005)
for i in pix_array_test1:
    for j in i :
        if float(j) < WindowCenter:
            pix_array_test1[k_i,k_j] = 0
        elif float(j) >= WindowCenter:
            pix_array_test1[k_i,k_j] = 1
        k_j += 1
    k_i += 1
    k_j = 0

F_Column = np.zeros(int(fileread_IM0005.Columns))
F_Rows = np.zeros(int(fileread_IM0005.Rows))
MID_Column = int(fileread_IM0005.Columns)/2
MID_Row = int(fileread_IM0005.Rows)/2
TESTWIDTH_Column = 800
TESTWIDTH_Rows = 800
CIRCLE_Diameter = 190.0         #The Diameter of ACR Phantom
R_square = (CIRCLE_Diameter)**2.0
X_square = (CIRCLE_Diameter-Scale_x_test1)**2.0
NARROW = int(2.0*((R_square - X_square)**0.5))      #The Least Length of circle chord if diameter is CIRCLE_DIAMETER
START_Column = int(MID_Column - (TESTWIDTH_Column/2))
START_Row = int(MID_Row - (TESTWIDTH_Rows/2))
f_i = 0
x = 0.
while START_Column < int(MID_Column + (TESTWIDTH_Column/2)) :
    for x in pix_array_test1[START_Column] :
        F_Column[f_i] += x
        f_i +=1
    f_i = 0
    START_Column += 1

F_Column_Non_Ghosting = np.flatnonzero(F_Column > NARROW)
Left_Pixel = F_Column_Non_Ghosting[0]
Right_Pixel = F_Column_Non_Ghosting[-1]
LR_Length = (int(Right_Pixel) - int(Left_Pixel))*Scale_x_test1
print('Actual Left-Right Length',LR_Length,'mm.')
Error_LR_test1 = abs(LR_Length-CIRCLE_Diameter)*100/CIRCLE_Diameter
print('Error = {:.4f} %'.format(Error_LR_test1))
if abs(LR_Length-190.0) > 2.0 :
    print('WARNING  :  Left-Right is not accuracy : Acceptable Value must be +- 2 mm from 190.0 mm.')
else :
    print('Left-Right Value is accuracy')

k_i = 0
for i in pix_array_test1:
    F_Rows[k_i] = np.sum(i)
    k_i += 1

F_Row_Non_Ghosting = np.flatnonzero(F_Rows > NARROW)
Top_Pixel = F_Row_Non_Ghosting[0]
Bottom_Pixel = F_Row_Non_Ghosting[-1]
TB_Length = (int(Bottom_Pixel) - int(Top_Pixel))*Scale_y_test1
print('Actual Top-Bottom Length',TB_Length,'mm.')
Error_TB_test1 = abs(TB_Length-CIRCLE_Diameter)*100/CIRCLE_Diameter
print('Error = {:.4f} %'.format(Error_TB_test1))
if abs(TB_Length-190.0) > 2.0 :
    print('WARNING  :  Top-Bottom is not accuracy : Acceptable Value must be +- 2 mm from 190.0 mm.')
else :
    print('Top-Bottom Value is accuracy')

Left_standard = 97
Right_standard = 704
Top_standard = 106
Bottom_standard = 710
Left_calibrating = int(Left_standard-Left_Pixel)
Right_calibrating = int(Right_standard-Right_Pixel)
Top_calibrating = int(Top_standard-Top_Pixel)
Bottom_calibrating = int(Bottom_standard-Bottom_Pixel)

#### Processing High Contrast Spatial Resolution Test
print()
print('*********  MRI QUALITY CONTROLLING : High Contrast Spatial Resolution Test  *********')
print()
print("Testing from DICOM\IM_0001")
Scale_x_test2 = fileread_IM0001.PixelSpacing[0]
Scale_y_test2 = fileread_IM0001.PixelSpacing[1]
Real_x_test2 = int(fileread_IM0001.Rows)*float(Scale_x_test2)
Real_y_test2 = int(fileread_IM0001.Columns)*float(Scale_y_test2)
print('IM_0001 Resolution = {0:} x {1:}'.format(fileread_IM0001.Rows,fileread_IM0001.Columns))
print('IM_0001 Image size = ',Real_x_test2,'mm','x',Real_y_test2,'mm')
WindowWidth_test2 = float(fileread_IM0001.WindowWidth)
WindowCenter_test2 = float(fileread_IM0001.WindowCenter)
print('IM_0001 Window Width = ',WindowWidth_test2)
print('IM_0001 Window Center = ',WindowCenter_test2)
print()

pix_array_test2 = np.copy(pix_array_IM_0001)
array_test2 = pix_array_test2[495-(Top_calibrating):560-(Bottom_calibrating),320-(Left_calibrating):530-(Right_calibrating)]
fig_test2 , ax_test2 = plt.subplots()
ax_test2.imshow(array_test2, cmap='gray')
ax_test2.axis('off')
ax_test2.set_title('Holes Arrays 1.1 1.0 and 0.9')
plt.show()

#### Processing Slice Thickness Accuracy Test
print()
print('*********  MRI QUALITY CONTROLLING : Slice Thickness Accuracy Test  *********')
print()
print("Testing from DICOM\IM_0001")
Scale_x_test3 = fileread_IM0001.PixelSpacing[0]
Scale_y_test3 = fileread_IM0001.PixelSpacing[1]
Real_x_test3 = int(fileread_IM0001.Rows)*float(Scale_x_test3)
Real_y_test3 = int(fileread_IM0001.Columns)*float(Scale_y_test3)
print('IM_0001 Resolution = {0:} x {1:}'.format(fileread_IM0001.Rows,fileread_IM0001.Columns))
print('IM_0001 Image size = ',Real_x_test3,'mm','x',Real_y_test3,'mm')
WindowWidth_test3 = float(fileread_IM0001.WindowWidth)
WindowCenter_test3 = float(fileread_IM0001.WindowCenter)
print('IM_0001 Window Width = ',WindowWidth_test3)
print('IM_0001 Window Center = ',WindowCenter_test3)
print()
pix_array_test3 = np.copy(pix_array_IM_0001)
array_ROI_T = pix_array_test3[390-(Top_calibrating):405-(Bottom_calibrating),150-(Left_calibrating):650-(Right_calibrating)]
array_ROI_B = pix_array_test3[405-(Top_calibrating):420-(Bottom_calibrating),150-(Left_calibrating):650-(Right_calibrating)]
Sensitivity = 13        #Adjust Level of Sensitivity from 1 - 15
print("Least Pixel Sensitivity is {} from 15 Pixel".format(Sensitivity))
k_i = 0
k_j = 0
MEAN_T = int(np.mean(array_ROI_T))
MEAN_B = int(np.mean(array_ROI_B))
HALF_T = MEAN_T/2
HALF_B = MEAN_B/2

for i in array_ROI_T:
    for j in i :
        if float(j) < HALF_T:
            array_ROI_T[k_i,k_j] = 0
        elif float(j) >= HALF_T:
            array_ROI_T[k_i,k_j] = 1
        k_j += 1
    k_i += 1
    k_j = 0

k_i = 0
k_j = 0
for i in array_ROI_B:
    for j in i :
        if float(j) < HALF_B:
            array_ROI_B[k_i,k_j] = 0
        elif float(j) >= HALF_B:
            array_ROI_B[k_i,k_j] = 1
        k_j += 1
    k_i += 1
    k_j = 0

S_T = np.zeros(500)
S_B = np.zeros(500)
k_j = 0
for i in array_ROI_T:
    for j in i :
        S_T[k_j] += float(j)
        k_j += 1
    k_j = 0
k_j = 0
for i in array_ROI_B:
    for j in i :
        S_B[k_j] += float(j)
        k_j += 1
    k_j = 0

S_T_L = np.flatnonzero(S_T > Sensitivity)
S_B_L = np.flatnonzero(S_B > Sensitivity)
Left_T = int(S_T_L[0])
Right_T = int(S_T_L[-1])
Top_Length = (Right_T - Left_T)*Scale_x_test3
Left_B = int(S_B_L[0])
Right_B = int(S_B_L[-1])
Bottom_Length = (Right_B - Left_B)*Scale_x_test3
print("Top Slice Length = {:.2f} mm".format(Top_Length))
print("Bottom Slice Length = {:.2f} mm".format(Bottom_Length))
Slice_Acc = (0.2*(Top_Length*Bottom_Length))/(Top_Length+Bottom_Length)
Error_test3 = abs(Slice_Acc-5.0)*100/5.0
print("Slice Thickness = {:.2f} mm".format(Slice_Acc))
print("Error = {:.4f} %".format(Error_test3))
if abs(Slice_Acc - 5.0) > 0.7:
    print("Warning: Slice Thickness is not accuracy, Value Must be +- 0.7 from 5.0 mm")

#### Processing Slice Position Accuracy Test
print()
print('*********  MRI QUALITY CONTROLLING : Slice Position Accuracy Test  *********')
print()
print("Testing from DICOM\IM_0001 and DICOM\IM_0011")
Scale_x1_test4 = fileread_IM0001.PixelSpacing[0]
Scale_y1_test4 = fileread_IM0001.PixelSpacing[1]
Real_x1_test4 = int(fileread_IM0001.Rows)*float(Scale_x1_test4)
Real_y1_test4 = int(fileread_IM0001.Columns)*float(Scale_y1_test4)
Scale_x2_test4 = fileread_IM0011.PixelSpacing[0]
Scale_y2_test4 = fileread_IM0011.PixelSpacing[1]
Real_x2_test4 = int(fileread_IM0011.Rows)*float(Scale_x2_test4)
Real_y2_test4 = int(fileread_IM0011.Columns)*float(Scale_y2_test4)
print('IM_0001 Resolution = {0:} x {1:}'.format(fileread_IM0001.Rows,fileread_IM0001.Columns))
print('IM_0001 Image Size = ',Real_x1_test4,'mm','x',Real_y1_test4,'mm')
print('IM_0011 Resolution = {0:} x {1:}'.format(fileread_IM0001.Rows,fileread_IM0001.Columns))
print('IM_0011 Image Size = ',Real_x2_test4,'mm','x',Real_y2_test4,'mm')
WindowWidth1_test4 = float(fileread_IM0001.WindowWidth)
WindowCenter1_test4 = float(fileread_IM0001.WindowCenter)
WindowWidth2_test4 = float(fileread_IM0011.WindowWidth)
WindowCenter2_test4 = float(fileread_IM0011.WindowCenter)
print('IM_0001 Window Width = ',WindowWidth1_test4)
print('IM_0001 Window Center = ',WindowCenter1_test4)
print('IM_0011 Window Width = ',WindowWidth2_test4)
print('IM_0011 Window Center = ',WindowCenter2_test4)

pix_array_IM_0001_test4 = np.copy(pix_array_IM_0001)
pix_array_IM_0011_test4 = np.copy(pix_array_IM_0011)
array_IM0001 = pix_array_IM_0001_test4[170-(Top_calibrating):220-(Bottom_calibrating),365-(Left_calibrating):435-(Right_calibrating)]
array_IM0011 = pix_array_IM_0011_test4[170-(Top_calibrating):220-(Bottom_calibrating),365-(Left_calibrating):435-(Right_calibrating)]

k_i = 0
k_j = 0
Signal_IM0001 = int(WindowCenter1_test4)/5
Signal_IM0011 = int(WindowCenter2_test4)/5
for i in array_IM0001 :
    for j in i :
        if float(j) < Signal_IM0001:
            array_IM0001[k_i,k_j] = 1
        elif float(j) >= Signal_IM0001:
            array_IM0001[k_i,k_j] = 0
        k_j += 1
    k_i += 1
    k_j = 0
k_i = 0
k_j = 0
for i in array_IM0011 :
    for j in i :
        if float(j) < Signal_IM0011:
            array_IM0011[k_i,k_j] = 1
        elif float(j) >= Signal_IM0011:
            array_IM0011[k_i,k_j] = 0
        k_j += 1
    k_i += 1
    k_j = 0

Value_IM0001_x = np.zeros(70)
Value_IM0011_x = np.zeros(70)
Value_IM0001_y = np.zeros(50)
Value_IM0011_y = np.zeros(50)
k_j = 0
k_i = 0
for i in array_IM0001:
    Value_IM0001_y[k_i] = np.sum(i)
    for j in i :
        Value_IM0001_x[k_j] += float(j)
        k_j += 1
    k_j = 0
    k_i += 1
k_j = 0
k_i = 0
for i in array_IM0011:
    Value_IM0011_y[k_i] = np.sum(i)
    for j in i :
        Value_IM0011_x[k_j] += float(j)
        k_j += 1
    k_j = 0
    k_i += 1
k_j = 0
k_i = 0
MAX_IM0001_x = np.max(Value_IM0001_x) -3        #Actual Value might be over than the real
MAX_IM0011_x = np.max(Value_IM0011_x) -3        #Actual Value might be over than the real
MAX_IM0001_y = np.max(Value_IM0001_y) -3
MAX_IM0011_y = np.max(Value_IM0011_y) -3
Value_IM0001_NonZero1 = np.flatnonzero(Value_IM0001_x>MAX_IM0001_x)
Value_IM0011_NonZero1 = np.flatnonzero(Value_IM0011_x>MAX_IM0011_x)
Value_IM0001_NonZero2 = np.flatnonzero(Value_IM0001_x>10)
Value_IM0011_NonZero2 = np.flatnonzero(Value_IM0011_x>10)
Value_IM0001_NonZero1_Y = np.flatnonzero(Value_IM0001_y>MAX_IM0001_y)
Value_IM0011_NonZero1_Y = np.flatnonzero(Value_IM0011_y>MAX_IM0011_y)
Value_IM0001_NonZero2_Y = np.flatnonzero(Value_IM0001_y)
Value_IM0011_NonZero2_Y = np.flatnonzero(Value_IM0011_y)
L_x_IM0001 = int(Value_IM0001_NonZero2[0])-int(Value_IM0001_NonZero1[0])-int(Value_IM0001_NonZero2[-1])+int(Value_IM0001_NonZero1[-1])
L_x_IM0011 = int(Value_IM0011_NonZero2[0])-int(Value_IM0011_NonZero1[0])-int(Value_IM0011_NonZero2[-1])+int(Value_IM0011_NonZero1[-1])
L_y_IM0001 = int(Value_IM0001_NonZero2_Y[-1])-int(Value_IM0001_NonZero1_Y[-1])
if L_x_IM0001 < 0.00 :
    L_y_IM0001 = (0 - L_y_IM0001)
L_y_IM0011 = int(Value_IM0011_NonZero2_Y[-1])-int(Value_IM0011_NonZero1_Y[-1])
if L_x_IM0011 < 0.00 :
    L_y_IM0011 = (0 - L_y_IM0011)
L_IM0001 = L_y_IM0001*Scale_x1_test4
L_IM0011 = L_y_IM0011*Scale_x2_test4
print("L from IM_0001 = {:.2f} mm".format(L_IM0001))
print("L from IM_0011 = {:.2f} mm".format(L_IM0011))
print("Slice Position Accuracy = {}".format(L_IM0011+L_IM0001))
if abs(L_IM0011+L_IM0001) > 5:
    print("WARNING : There are some Slice Position inaccuracy")

#### Processing Percent Image Uniformity
print()
print('*********  MRI QUALITY CONTROLLING : Percent Image Uniformity Test  *********')
print("Testing from DICOM\IM_0007")
Scale_x_test5 = fileread_IM0007.PixelSpacing[0]
Scale_y_test5 = fileread_IM0007.PixelSpacing[1]
Real_x_test5 = int(fileread_IM0007.Rows)*float(Scale_x_test5)
Real_y_test5 = int(fileread_IM0007.Columns)*float(Scale_y_test5)
print('IM_0007 Resolution = {0:} x {1:}'.format(fileread_IM0005.Rows,fileread_IM0005.Columns))
print('IM_0007 Image size = ',Real_x_test5,'mm','x',Real_y_test5,'mm')
WindowWidth_test5 = float(fileread_IM0007.WindowWidth)
WindowCenter_test5 = float(fileread_IM0007.WindowCenter)
print('IM_0007 Window Width = ',WindowWidth_test5)
print('IM_0007 Window Center = ',WindowCenter_test5)
print()

pix_array_test5 = np.copy(pix_array_IM_0007)
WindowCenter_test5 = WindowWidth_test5/2

def Window_Adjust(array,windowwidth1,windowwidth2):
    k_i = 0
    k_j = 0
    new_array = np.zeros_like(array)
    for i in array:
        for j in i :
            new_array[k_i,k_j] = j*(windowwidth2/windowwidth1)
            k_j += 1
        k_j = 0
        k_i += 1
    return new_array

def Window_Uniform(array,windowwidth1,windowlevel,windowwidth2):
    k_i = 0
    k_j = 0
    new_array = np.zeros_like(array)
    k = 0
    max_k = windowlevel + (windowwidth1/2)
    min_k = windowlevel - (windowwidth1/2)
    for i in array:
        for j in i :
            if float(j) < min_k :
                new_array[k_i,k_j] = 0
            if float(j) > max_k :
                new_array[k_i,k_j] = windowwidth2
            if float(j) >= min_k and float(j) <= max_k :
                new_array[k_i,k_j] = (j-min_k)*(windowwidth2/windowwidth1)
            k_j += 1
        k_j = 0
        k_i += 1
    return new_array

WindowWidth2_test5 = 256
WindowLevel_test5 = WindowCenter_test5
cal = 0
print("Finding Low Signal Value...")
while cal < 50000 or cal > 200000:
    cal = 0
    test_array_test5 = Window_Uniform(pix_array_test5,1,WindowLevel_test5,WindowWidth2_test5)
    k_i = 0
    k_j = 0
    for i in test_array_test5:
        for j in i:
            if float(j) > 0:
                cal += 1
            k_j += 1
        k_j = 0
        k_i += 1
    print("Current White Area {:.2f} (mm^2)".format(cal*Scale_x_test5*Scale_y_test5))
    if cal > 210000:
        WindowLevel_test5 += 10
    if cal < 210000:
        WindowLevel_test5 += 1
signal_low = WindowLevel_test5
test_array2_test5 = Window_Uniform(pix_array_test5,1,WindowLevel_test5,WindowWidth2_test5)
print()
print("Finding High Signal Value...")
while cal > 1024:
    cal = 0
    test_array2_test5 = Window_Uniform(pix_array_test5,1,WindowLevel_test5,WindowWidth2_test5)
    k_i = 0
    k_j = 0
    for i in test_array2_test5:
        for j in i:
            if float(j) > 0:
                cal += 1
            k_j += 1
        k_j = 0
        k_i += 1
    print("Current White Area {:.2f} (mm^2)".format(cal*Scale_x_test5*Scale_y_test5))
    if cal > 7000:
        WindowLevel_test5 += 10
    if cal < 7000:
        WindowLevel_test5 += 1

signal_high = WindowLevel_test5
print("")
print("Low Signal = {}".format(signal_low))
print("High Signal = {}".format(signal_high))
percent_uniform = (1-((signal_high-signal_low)/(signal_low+signal_high)))*100
print("PIU = {:.4f} %".format(percent_uniform))
tesla = float(fileread_IM0007.MagneticFieldStrength)
print("MRI Magnetic Field Strength = {} T".format(tesla))
if percent_uniform < 82.0 and tesla <= 3.0 :
    print("Warning : Low PIU")
if percent_uniform < 87.5 and tesla > 3.0 :
    print("Warning : Low PIU")

#### Processing Low Contrast Detectability Test
print()
print('*********  MRI QUALITY CONTROLLING : Low Contrast Detectability Test  *********')
print("Testing from DICOM\IM_0008 , IM_0009 , IM_0010 and IM_0011")
Scale_x_IM0008_test6 = fileread_IM0008.PixelSpacing[0]
Scale_y_IM0008_test6 = fileread_IM0008.PixelSpacing[1]
Scale_x_IM0009_test6 = fileread_IM0009.PixelSpacing[0]
Scale_y_IM0009_test6 = fileread_IM0009.PixelSpacing[1]
Scale_x_IM0010_test6 = fileread_IM0010.PixelSpacing[0]
Scale_y_IM0010_test6 = fileread_IM0010.PixelSpacing[1]
Scale_x_IM0011_test6 = fileread_IM0011.PixelSpacing[0]
Scale_y_IM0011_test6 = fileread_IM0011.PixelSpacing[1]
Real_x_IM0008_test6 = int(fileread_IM0008.Rows)*float(Scale_x_IM0008_test6)
Real_y_IM0008_test6 = int(fileread_IM0008.Columns)*float(Scale_y_IM0008_test6)
Real_x_IM0009_test6 = int(fileread_IM0009.Rows)*float(Scale_x_IM0009_test6)
Real_y_IM0009_test6 = int(fileread_IM0009.Columns)*float(Scale_y_IM0009_test6)
Real_x_IM0010_test6 = int(fileread_IM0010.Rows)*float(Scale_x_IM0010_test6)
Real_y_IM0010_test6 = int(fileread_IM0010.Columns)*float(Scale_y_IM0010_test6)
Real_x_IM0011_test6 = int(fileread_IM0011.Rows)*float(Scale_x_IM0011_test6)
Real_y_IM0011_test6 = int(fileread_IM0011.Columns)*float(Scale_y_IM0011_test6)
print('IM_0008 Resolution = {0:} x {1:}'.format(fileread_IM0008.Rows,fileread_IM0008.Columns))
print('IM_0008 Image size = ',Real_x_IM0008_test6,'mm','x',Real_y_IM0008_test6,'mm')
print('IM_0009 Resolution = {0:} x {1:}'.format(fileread_IM0009.Rows,fileread_IM0009.Columns))
print('IM_0009 Image size = ',Real_x_IM0009_test6,'mm','x',Real_y_IM0009_test6,'mm')
print('IM_0010 Resolution = {0:} x {1:}'.format(fileread_IM0010.Rows,fileread_IM0010.Columns))
print('IM_0010 Image size = ',Real_x_IM0010_test6,'mm','x',Real_y_IM0010_test6,'mm')
print('IM_0011 Resolution = {0:} x {1:}'.format(fileread_IM0011.Rows,fileread_IM0011.Columns))
print('IM_0011 Image size = ',Real_x_IM0011_test6,'mm','x',Real_y_IM0011_test6,'mm')
WindowWidth_IM0008_test6 = float(fileread_IM0008.WindowWidth)
WindowCenter_IM0008_test6 = float(fileread_IM0008.WindowCenter)
WindowWidth_IM0009_test6 = float(fileread_IM0009.WindowWidth)
WindowCenter_IM0009_test6 = float(fileread_IM0009.WindowCenter)
WindowWidth_IM0010_test6 = float(fileread_IM0010.WindowWidth)
WindowCenter_IM0010_test6 = float(fileread_IM0010.WindowCenter)
WindowWidth_IM0011_test6 = float(fileread_IM0011.WindowWidth)
WindowCenter_IM0011_test6 = float(fileread_IM0011.WindowCenter)
print('IM_0008 Window Width = ',WindowWidth_IM0008_test6)
print('IM_0008 Window Center = ',WindowCenter_IM0008_test6)
print('IM_0009 Window Width = ',WindowWidth_IM0009_test6)
print('IM_0009 Window Center = ',WindowCenter_IM0009_test6)
print('IM_0010 Window Width = ',WindowWidth_IM0010_test6)
print('IM_0010 Window Center = ',WindowCenter_IM0010_test6)
print('IM_0011 Window Width = ',WindowWidth_IM0011_test6)
print('IM_0011 Window Center = ',WindowCenter_IM0011_test6)
print()
pix_array_IM0008_test6 = np.copy(pix_array_IM_0008)
pix_array_IM0009_test6 = np.copy(pix_array_IM_0009)
pix_array_IM0010_test6 = np.copy(pix_array_IM_0010)
pix_array_IM0011_test6 = np.copy(pix_array_IM_0011)
test_array_IM0008_test6 = pix_array_IM0008_test6[279-(Top_calibrating):581-(Bottom_calibrating),259-(Left_calibrating):561-(Right_calibrating)]
test_array_IM0009_test6 = pix_array_IM0009_test6[279-(Top_calibrating):581-(Bottom_calibrating),259-(Left_calibrating):561-(Right_calibrating)]
test_array_IM0010_test6 = pix_array_IM0010_test6[279-(Top_calibrating):581-(Bottom_calibrating),259-(Left_calibrating):561-(Right_calibrating)]
test_array_IM0011_test6 = pix_array_IM0011_test6[279-(Top_calibrating):581-(Bottom_calibrating),259-(Left_calibrating):561-(Right_calibrating)]
parameter1_IM0008 = int(WindowWidth_IM0008_test6*0.2101)
parameter2_IM0008 = int(WindowCenter_IM0008_test6*1.3052)
test_array_IM0008_test6 = Window_Uniform(test_array_IM0008_test6,parameter1_IM0008,parameter2_IM0008,WindowWidth_IM0008_test6)
parameter1_IM0009 = int(WindowWidth_IM0009_test6*0.2101)
parameter2_IM0009 = int(WindowCenter_IM0009_test6*1.3052)
test_array_IM0009_test6 = Window_Uniform(test_array_IM0009_test6,parameter1_IM0009,parameter2_IM0009,WindowWidth_IM0009_test6)
parameter1_IM0010 = int(WindowWidth_IM0010_test6*0.2101)
parameter2_IM0010 = int(WindowCenter_IM0010_test6*1.3052)
test_array_IM0010_test6 = Window_Uniform(test_array_IM0010_test6,parameter1_IM0010,parameter2_IM0010,WindowWidth_IM0010_test6)
parameter1_IM0011 = int(WindowWidth_IM0011_test6*0.2101)
parameter2_IM0011 = int(WindowCenter_IM0011_test6*1.3052)
test_array_IM0011_test6 = Window_Uniform(test_array_IM0011_test6,parameter1_IM0011,parameter2_IM0011,WindowWidth_IM0011_test6)
fig1, ax1 = plt.subplots()
ax1.imshow(test_array_IM0008_test6, cmap="gray")
ax1.axis('off')
fig1.suptitle("IM_0008")
fig2, ax2 = plt.subplots()
ax2.imshow(test_array_IM0009_test6, cmap="gray")
ax2.axis('off')
fig2.suptitle("IM_0009")
fig3, ax3 = plt.subplots()
ax3.imshow(test_array_IM0010_test6, cmap="gray")
ax3.axis('off')
fig3.suptitle("IM_0010")
fig4, ax4 = plt.subplots()
ax4.imshow(test_array_IM0011_test6, cmap="gray")
ax4.axis('off')
fig4.suptitle("IM_0011")
plt.show()

#### Processing Percent Signal Ghosting Test
print()
print('********  MRI QUALITY CONTROLLING : Percent Signal Ghosting Test  *********')

pix_array_test7_8 = np.copy(pix_array_IM_0007)
ROI = pix_array_test7_8[190-(Top_calibrating):590-(Bottom_calibrating),190-(Left_calibrating):590-(Right_calibrating)]
mean_roi = np.mean(ROI)
ROI_01 = pix_array_test7_8[350-(Top_calibrating):450-(Bottom_calibrating),350-(Left_calibrating):450-(Right_calibrating)]
mean_roi01 = np.mean(ROI_01)
ROI_02 = pix_array_test7_8[200-(Top_calibrating):300-(Bottom_calibrating),500-(Left_calibrating):600-(Right_calibrating)]
mean_roi02 = np.mean(ROI_02)
Top_ROI = pix_array_test7_8[10-(Top_calibrating):55-(Bottom_calibrating),150-(Top_calibrating):650-(Right_calibrating)]
mean_top = np.mean(Top_ROI)
std_top = np.std(Top_ROI)
Left_ROI = pix_array_test7_8[150-(Top_calibrating):650-(Bottom_calibrating),10-(Left_calibrating):55-Right_calibrating]
mean_left = np.mean(Left_ROI)
std_left = np.std(Left_ROI)
Right_ROI = pix_array_test7_8[150-(Top_calibrating):650-(Bottom_calibrating),745-(Left_calibrating):790-(Right_calibrating)]
mean_right = np.mean(Right_ROI)
std_right = np.std(Right_ROI)
Bottom_ROI = pix_array_test7_8[745-(Top_calibrating):790-(Bottom_calibrating),150-(Left_calibrating):650-(Right_calibrating)]
mean_bottom = np.mean(Bottom_ROI)
std_bottom = np.std(Bottom_ROI)

print("Large ROI Mean Signal = {:.2f}".format(mean_roi))
print("Top ROI Mean Signal = {:.2f}".format(mean_top))
print("Left ROI Mean Signal = {:.2f}".format(mean_left))
print("Right ROI Mean Signal = {:.2f}".format(mean_right))
print("Bottom ROI Mean Signal = {:.2f}".format(mean_bottom))

ghosting_ratio = ((mean_top+mean_bottom)-(mean_left+mean_right))/(2*mean_roi)

print("Percent Signal Ghosting Ratio = {:.5f}".format(ghosting_ratio))

if ghosting_ratio > 0.025 :
    print("Warning : High Percent Signal Ghosting, value must below 0.0025 ")

#### Processing Signal to Noise Ratio Test
print()
print('********  MRI QUALITY CONTROLLING : Signal to Noise Test  *********')

Mean_X = (mean_roi01+mean_roi02)/2
Mean_Y = (mean_bottom+mean_right+mean_left+mean_top)/4
Mean_SD = (std_bottom+std_right+std_left+std_top)/4
SNR = (Mean_X-Mean_Y)/Mean_SD
print("Center ROI Mean Signal = {:.2f}".format(mean_roi01))
print("Border ROI Mean Signal = {:.2f}".format(mean_roi02))
print("Top ROI Mean Signal = {:.2f}".format(mean_top))
print("S.D. of Top Background = {:.2f}".format(std_top))
print("Left ROI Mean Signal = {:.2f}".format(mean_left))
print("S.D. of Left Background = {:.2f}".format(std_left))
print("Right ROI Mean Signal = {:.2f}".format(mean_right))
print("S.D. of Right Background = {:.2f}".format(std_right))
print("Bottom ROI Mean Signal = {:.2f}".format(mean_bottom))
print("S.D. of Bottom Background = {:.2f}".format(std_bottom))

print("Signal to Noise Ratio = {:.4f}".format(SNR))
end_time = time.time()
print("**************************************END*************************************")
print("Total Execution Time = {:.4f} seconds".format(end_time-start_time))
