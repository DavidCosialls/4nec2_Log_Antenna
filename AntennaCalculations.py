import math

import matplotlib.pyplot as plt
import numpy as np 
from tkinter import *
from tkinter import messagebox

from PIL import ImageTk, Image

f = open('logAntenna', 'w')

f.write("CM\n")
f.write("CE\n")

def openImpedanceWindow(): 
    newWindow = Toplevel(master)

    newWindow.title("Graph")

    # sets the geometry of toplevel
    newWindow.geometry("400x400")

    # A Label widget to show in toplevel
    # Create an object of tkinter ImageTk

    img = ImageTk.PhotoImage(Image.open("impedance.jpg"))
    # Create a Label Widget to display the text or Image

    label= Label(newWindow, image= img)
    label.image= img
    label.pack()



master = Tk()
master.title("SRF log Antenna Calculations")

master.geometry('1050x600')
#DEFINITION OF THE BUTTONS AND LABELS
titleLbl = Label(master, text="Select SIGMA and THAO for the desired Directivity")
titleLbl.place(x=0, y=0)

maxfLbl = Label(master, text = "Maximum frequency (MHz): ")
maxfLbl.place(x=0, y=40)
maxfLbl_entry = Entry(master, width=10)
maxfLbl_entry.place(x=180, y=40)

minfLbl = Label(master, text = "Minimum frequency (MHz): ")
minfLbl.place(x=0, y=80)
minfLbl_entry = Entry(master, width=10)
minfLbl_entry.place(x=180, y=80)

sigmaLbl = Label(master, text = "Sigma Value from the graph: ")
sigmaLbl.place(x=0, y=120)
sigma_entry = Entry(master,width=10)
sigma_entry.place(x=180, y=120)
#sigma = float(sigma_entry.get())

thaoLbl = Label(master, text = "Thao value from the graph: ")
thaoLbl.place(x=0, y=160)
thao_entry = Entry(master, width=10)
thao_entry.place(x=180, y=160)

widthLbl = Label(master, text = "Wire width value (m): ")
widthLbl.place(x=0, y=200)
width_entry = Entry(master, width=10)
width_entry.place(x=180, y=200)

RinLbl = Label(master, text = "Rin (Ohms): ")
Rin_entry = Entry(master, width=10)
RinLbl.place(x=0, y=240)
Rin_entry.place(x=180, y=240)

SCoeffLbl = Label(master, text = "Spacing Coeff (default = 1)")
SCoeff_entry = Entry(master, width=10)
SCoeffLbl.place(x=0, y=280)
SCoeff_entry.place(x=180, y=280)


img = ImageTk.PhotoImage(Image.open("directivityGRAPH.jpg"))

directivityGRAPHLBL = Label(master, image= img)
directivityGRAPHLBL.image= img
directivityGRAPHLBL.place(x=330, y=10)

#FINISH OF THE DEFINITION

c = 3.0*10**8 #m/s
fmax = 1600.0*10**6 #Hz
fmin = 1200.0*10**6 #Hz

alpha = 0.0 
B = 0.0
Bar = 0.0
Bs = 0.0
lambdaMax = 0.0
L = 0.0
N = 0
lmax = 0.0
lmaxOverdmax = 0.0
sigma = 0.0
thao = 0.0
dmax = 0.0
Rin = 0.0
Za = 0.0
ZaRin = 0.0
sCoeff = 0.0

def calculations():
    
    DesiredDirectivity = 10.5 #dB

    done = True
    
    #GLOBAL VARIABLES
    global alpha
    global B 
    global Bar 
    global Bs 
    global lambdaMax
    global L
    global N 
    global lmax 
    global lmaxOverdmax 
    global sigma
    global thao 
    global dmax 
    global Rin 
    global Za 
    global ZaRin 
    global dmax
    global sCoeff

    #TRY EVERY VARIABLE TO CHECK IF IT IS A CORRECT ENTRY VALUE

    try:
        fmax = float(maxfLbl_entry.get())*10**6
        fmin = float(minfLbl_entry.get())*10**6
    except:
        messagebox.showerror(title="ERROR",message="Sigma and Thao are not correctly introduced")
        return

    try:
        thao = float(thao_entry.get())
        sigma = float(sigma_entry.get())
    except:
        messagebox.showerror(title="ERROR",message="Sigma and Thao are not correctly introduced")
        return 

    try: 
        dmax = float(width_entry.get())
    except:
        messagebox.showerror(title="ERROR", message="The width is not correctly introduced")
        return

    if(thao < 0 or thao >= 1.0):
        messagebox.showerror(title="ERROR",message="Incorrect value for Thao")
        return 
    
    if (sigma > 0.22 or sigma<0.04):
        messagebox.showerror(title="ERROR",message="Incorrect value for Sigma")
        return 

    try: 
        Rin = float(Rin_entry.get())
    except:
        messagebox.showerror(title="ERROR",message="Incorrect value for Rin")

    try:
        sCoeff = float(SCoeff_entry.get())
    except:
        messagebox.showerror(title="ERROR", message="Incorrect value for S_Coefficient, please value from 0.1 to 100")

    #FINISH OF DEFINING VARIABLES

    sigmaOpt = 0.243 * thao - 0.051 #WE CHOOSE THE OPTIMUM SIGMA
    sigma = sigmaOpt
    sigma_entry.config(textvariable=str(sigma))

    alpha = np.arctan( (1.0-thao)/(4.0*sigma) ) #HOW TO CALCULATE ALPHA

    cotAlpha = 1/(np.tan(alpha)) #COTANGENT

    Bar = 1.1+7.7*(1.0-thao)**2 *(1.0/(np.tan(alpha))) #BANDWIDTH
    
    B = fmax/fmin 
    Bs = B*Bar

    lambdaMax = c/fmin

    L = (lambdaMax/4.0)*(1.0-(1/Bs))*(1.0/np.tan(alpha))

    N = (1.0 + (np.log(Bs)/np.log(1.0/thao))) #NUMBER OF ELEMENTS OF THE ANTENNA.
    N = np.round(N)

    lmax = lambdaMax/2.0  

    ell_Z_term = 0.125 * lambdaMax #Z TERM TO
    ell_Z_term = np.round(ell_Z_term, 4)

    feed_length = 0.1

    #WRITE ALL THE VARIABLES ON THE FILE
    f.write("SY\tfmin="+str(int(fmin))+"\n")
    f.write("SY\tfmax="+str(int(fmax))+"\n")
    f.write("SY\talpha="+str(np.rad2deg(alpha))+"\n")
    f.write("SY\tcotalpha="+str(cotAlpha)+"\n")
    f.write("SY\tc="+str(3*10**8)+"\n")
    f.write("SY\tB="+"fmax/fmin"+"\n")
    f.write("SY\tthao="+str(thao)+"\n")
    f.write("SY\tsigma="+str(sigma)+"\n")
    f.write("SY\tN="+str(int(N))+"\n")
    f.write("SY\tdmax="+str(dmax/2)+"\n")
    f.write("SY\tZ_term="+str(ell_Z_term)+"\n")
    f.write("SY\tfeed_length="+str(feed_length)+"\n")
    f.write("SY\tlambdaMax="+str(lambdaMax)+"\n")
    f.write("SY\tsegments=5"+"\n")

    Dipolelengths = []
    Dipolelengths_STRING = []
    i = 0
    while (i<N):
        Dipolelengths.append((lambdaMax*0.5)*(thao)**i)
        Dipolelengths_STRING.append("lambdaMax*0.5*thao^"+str(i))
        i += 1

    i=0
    d = []
    BoomLength = 0
    Xcoor = []
    d_String = []
    Xcoor_String = []
    while (i<N-1):
        d.append(0.5 * (Dipolelengths[i] - Dipolelengths[i+1]) * cotAlpha)
        d_String.append("0.5 *("+Dipolelengths_STRING[i]+"-"+Dipolelengths_STRING[i+1]+")*cotalpha")
        BoomLength += d[i]
        j = 0
        finalString = ""
        while (j <= i):
            if (j<i):
                finalString = finalString + d_String[j] + "+"
            if (j==i):
                finalString = finalString + d_String[j]
            j += 1

        Xcoor_String.append(finalString)
        Xcoor.append(BoomLength)
        i += 1       

    lmaxOverdmax = (Dipolelengths[int(N-1)])/dmax

    Za = 120.0*(np.log(lmaxOverdmax)-2.25)
    sigmaPrima = sigma/np.sqrt(thao)
    
    ZaRin = Za/Rin

    sigma_mean = sigma / np.sqrt(thao)
    Zc_feed = Rin**2 * (8 * sigma_mean * Za)**-1
    Zc_feed += Rin * np.sqrt( (Rin**2 * 0.015625 * sigma_mean**-2 * Za**-2) +1)

    

    print("\033[4;30;47m"+"---------RESULTS---------"+'\033[0;m')

    print("\033[1;33m"+"Alpha ="+'\033[0;m',np.rad2deg(alpha), "º")
    print("\033[1;33m"+"CotAlpha ="+'\033[0;m', cotAlpha)
    print("\033[1;33m"+"B ="+'\033[0;m',B)
    print("\033[1;33m"+"Bar ="+'\033[0;m',Bar)
    print("\033[1;33m"+"Bs ="+'\033[0;m',Bs)
    print("\033[1;33m"+"L ="+'\033[0;m',L)
    print("\033[1;33m"+"N ="+'\033[0;m',N)
    print("\033[1;33m"+"Sigma OPT ="+'\033[0;m',sigma)

    print("\033[1;33m"+"DIPOLE_LENGTHS : "+'\033[0;m')

    i=0
    s = dmax * np.cosh(Zc_feed/120)
    f.write("SY\tsCoeff=" + str(sCoeff) + "\n")
    f.write("SY\ts="+str(np.round(s, 4))+"\n")
    f.write("SY\tZTLINE="+str(np.round(Zc_feed, 4))+"\n")

    tag = 1
    segmentsVector = []
    while (i<len(Dipolelengths)):
        seg = np.floor(Dipolelengths[i] / s)
        if (seg % 2 == 0):
            seg = int(seg + 1)

        # MID POINT TO PUT THE T-LINE

        round = int(seg / 2) + 1
        segmentsVector.append(str(round))
        seg_string = str(int(seg))
        if (i==0):

            #                 TAG    SEGMENTS               X1                          Y1    Z1       X2                                           Y2                       Z2   RADIUS    COMMENT       
            #f.write("GW\t"+str(tag)+"\t"+str(int(seg))+"\t-feed_length\t-0.5*(lambdaMax*0.5)*(thao)^"+str(i)+"\t0\t-feed_length\t0.5*(lambdaMax*0.5)*(thao)^"+str(i)+"\t0\tdmax\t'DipoleUpper"+str(int(N-i))+"\n")

            f.write("GW\t"+str(tag)+"\t"+seg_string+"\t-feed_length\t0\t-0.5*(lambdaMax*0.5)*(thao)^"+str(i)+"\t-feed_length\t0\t0.5*(lambdaMax*0.5)*(thao)^"+str(i)+"\tdmax\t'Dipole number"+str(int(N-i))+"\n")

            tag += 1
            #f.write("GW\t"+str(tag)+"\tsegments\t-feed_length\t-s/2\t0\t-feed_length\t-0.5*(lambdaMax*0.5)*(thao)^"+str(i)+"\t0\tdmax\t'DipoleDown"+str(int(N-i))+"\n")
            #tag += 1
        else:
            if (i%2 != 0):
                #                 TAG    SEGMENTS               X1                          Y1    Z1       X2                                           Y2                       Z2   RADIUS    COMMENT
                f.write("GW\t" + str(tag) + "\t" + seg_string + "\t" +
                    Xcoor_String[i - 1] + "-feed_length\t0\t0.5*(lambdaMax*0.5)*(thao)^" + str(i) + "\t" +
                    Xcoor_String[i - 1] + "-feed_length\t0\t-0.5*(lambdaMax*0.5)*(thao)^" + str(i) + "\tdmax\t'Dipole number" + str(int(N - i)) + "\n")
                tag += 1
                # f.write("GW\t"+str(tag)+"\tsegments\t"+str(Xcoor[i-1])+"-feed_length\t-s/2\t0\t"+str(Xcoor[i-1])+"-feed_length\t-0.5*(lambdaMax*0.5)*(thao)^"+str(i)+"\t0\tdmax\t'DipoleDown"+str(int(N-i))+"\n")
                # tag += 1
            else:
                #                 TAG    SEGMENTS               X1                          Y1    Z1       X2                                           Y2                       Z2   RADIUS    COMMENT
                f.write("GW\t"+str(tag)+"\t"+seg_string+"\t" +
                        Xcoor_String[i-1]+"-feed_length\t0\t-0.5*(lambdaMax*0.5)*(thao)^"+str(i)+"\t" +
                        Xcoor_String[i-1]+"-feed_length\t0\t0.5*(lambdaMax*0.5)*(thao)^"+str(i)+"\tdmax\t'Dipole number "+str(int(N-i))+"\n")
                tag += 1
                #f.write("GW\t"+str(tag)+"\tsegments\t"+str(Xcoor[i-1])+"-feed_length\t-s/2\t0\t"+str(Xcoor[i-1])+"-feed_length\t-0.5*(lambdaMax*0.5)*(thao)^"+str(i)+"\t0\tdmax\t'DipoleDown"+str(int(N-i))+"\n")
                #tag += 1


        print("l"+str(int(i+1))+" =", Dipolelengths[i], "    s =", s, " segments = ", seg, "segments/2 = ", segmentsVector[i])
        i += 1

    f.write("GW\t"+str(tag)+"\t1\t-feed_length-Z_term\t0\ts/2\t-feed_length-Z_term\t0\t-s/2\tdmax\t'DownZ_term\n")
    tag += 1
    f.write("GW\t"+str(tag)+"\t1\t"+str(Xcoor[len(Xcoor)-1])+"\t0\t-2*dmax\t"+str(Xcoor[len(Xcoor)-1])+"\t0\t2*dmax\tdmax\t'FeedLine\n")
    print("")
    print("\033[1;33m"+"Distances: "+'\033[0;m')

    i = 0
    #SPACING
    Y_Spacing = []
   

    #      VOLTAGE SOURCE


    #        TAG 
    #EX	0	70	2	0	1	0	0
    f.write("GE\t0\n")
    f.write("GN\t-1\n")
    f.write("EK\n")
    #     Voltage TAG SEG OPT REAL IMAG MAGN PHASE 
    f.write("EX\t0\t"+str(tag)+"\t1\t0\t1\t0\t1\t0\t0\n")


    #Z-LINE
    f.write("TL\t"+str(tag)+"\t1\t"+str(tag-2)+"\t"+segmentsVector[tag-3]+"\tZTLINE\t0\t0\t0\t0\t0\n")
    #f.write("TL\t"+str(tag)+"\t3\t"+str(tag-3)+"\t1\t50\t0\t0\t0\t0\t0\n")

    #VOLTAGE SOURCE
    f.write("TL\t"+str(tag-1)+"\t1\t"+str(1)+"\t"+segmentsVector[0]+"\tZTLINE\t0\t0\t0\t0\t0\n")
    #f.write("TL\t"+str(tag-1)+"\t3\t"+str(2)+"\t1\t50\t0\t0\t0\t0\t0\n")

    tag = 1
    while (i<len(d)):
        
        Y_Spacing.append(s)
        print ("l"+str(int(i+1))+" to l"+str(int(i+2))+" = ", d[i], " @ Xcoor = ", Xcoor[i])


        #TL SEG_1 ANCHO SEG_2 ANCHO Z0 0 1e+99 1e+99 1e+99 1e+99
        #TL 1 s_w 2 s_w 50 0 1e+99 1e+99 1e+99 1e+99

        f.write("TL\t"+str(tag)+"\t"+segmentsVector[i]+"\t"+str(tag+1)+"\t"+segmentsVector[i+1]+"\tZTLINE\t0\t0\t0\t0\t0\n")
        i += 1
        tag += 1

    
    #print("\033[1;33m"+"Xcoor ="+'\033[0;m',Xcoor)
    print("\033[1;33m"+"Boom Length = "+'\033[0;m',BoomLength)
    print("\033[1;33m"+"Z_term = "+'\033[0;m', ell_Z_term)
    print("\033[1;33m"+"Za = "+'\033[0;m',Za)
    print("\033[1;33m"+"Zc feed = "+'\033[0;m',Zc_feed)
    print("\033[1;33m"+"s = "+'\033[0;m', s)

    AlphaText = "Alpha = "+str(np.rad2deg(alpha))
    AlphaLbl = Label(master, text = AlphaText)
    AlphaLbl.place(x=0, y=360)

    NText = "N = "+ str(N)
    NLbl = Label(master, text = NText)
    NLbl.place(x=0, y=400)

    
    f.write("FR\t0\t0\t0\t0\t1400\n")
    f.write("EN\n")
    f.close()


#
#lbl = Label(window, text="Hello")
#lbl.grid(column=0, row=0)
btn = Button(master, text="Calculate", command=calculations)
btn.place(x=95, y=320)









#MessageBox Example
#messagebox.showinfo('Message title', 'Message content')

master.mainloop()
    
