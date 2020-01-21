from ROOT import gROOT, TCanvas, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor
import ROOT as root
from math import cos, sin, pi, degrees, radians, pow, sqrt,acos
from itertools import permutations 
import sys
from array import array

def DeltaPhi( Phi1, Phi2 ):
    deltaphi=abs(Phi1-Phi2)
    if (deltaphi>pi):
        deltaphi=2*pi-deltaphi        
    return deltaphi

def DeltaPhiDir( Phi1, Phi2 ):
    deltaphi=Phi1-Phi2
    if(deltaphi>pi):
        deltaphi=deltaphi-2*pi   
    if(deltaphi<(-pi)):
        deltaphi=2*pi+deltaphi
    return deltaphi


def CalculatePerformance(pTH1F):
  #//expects an average of two quantities --> originally from total event energy, thus multiply by a factor of sqrt(2) later
    FLOAT_MAX=sys.float_info.max
    print 'float max',FLOAT_MAX
    resolution=0 
    resolutionError=0
    mean=0
    meanError=0
    
    
    if (5 > pTH1F.GetEntries()):    
        print pTH1F.GetName(),pTH1F.GetEntries()," entries) - skipped"
        return None;
    
    #Calculate raw properties of distribution
    sum = 0.
    total = 0.
    sx = 0
    sxx = 0.
    nbins=pTH1F.GetNbinsX()
    
    for i in range(0,nbins+1):      
        binx=pTH1F.GetBinLowEdge(i) + (0.5 * pTH1F.GetBinWidth(i));
        yi=pTH1F.GetBinContent(i)
        sx += yi * binx;
        sxx += yi * binx * binx;
        total += yi;
        
        
    rawMean=sx / total
    rawMeanSqua=sxx / total
    rawRms=sqrt(rawMeanSqua - rawMean * rawMean)
        
    sum0 = 0.
    is0 = 0
        
    
    for i in range(0,nbins+1) :
        if (sum0<(total / 10.)) :
            sum0 += pTH1F.GetBinContent(i);
            is0 = i;
           
  
    # Calculate truncated properties
    rmsmin=FLOAT_MAX
    sigma=FLOAT_MAX
    sigmasigma=FLOAT_MAX
    frac=FLOAT_MAX
    efrac=FLOAT_MAX
    mean=FLOAT_MAX
    low=FLOAT_MAX
    rms=FLOAT_MAX
    high=0.
  
    for istart in range(0,is0+1):
        sumn = 0.
        csum = 0.
        sumx = 0.
        sumxx = 0.
        iend = 0;
        
        for i in range (istart,nbins+1):
            if csum < (0.9 * total):
                binx=pTH1F.GetBinLowEdge(i) + (0.5 * pTH1F.GetBinWidth(i))
                yi=pTH1F.GetBinContent(i)
                csum += yi
                
                if (sumn < (0.9 * total)):
                    sumn += yi;
                    sumx += yi * binx;
                    sumxx+= yi * binx * binx;
                    iend = i;
                    
                   
      
        localMean=sumx / sumn
        localMeanSqua=sumxx / sumn
        localRms=sqrt(localMeanSqua - localMean * localMean)
	
        if (localRms < rmsmin):
	  
            mean = localMean;
            rms = localRms;
            low = pTH1F.GetBinLowEdge(istart);
            high = pTH1F.GetBinLowEdge(iend);
            rmsmin = localRms;
	    #mean=1.;
	    frac = rms / mean* 100.;
	    efrac = frac / sqrt(total);
	    meanError = mean / sqrt(total); 
	  
	resolution = frac;
	resolutionError = efrac;
        
    print "resolution/error ",resolution,resolutionError," mean ",mean,meanError
    return resolution,resolutionError,mean,meanError


def setUpperCanvas(canvas_name) :
    c1= TCanvas(canvas_name,canvas_name,10,50,600,500)
    c1.cd()
    gPad.SetTopMargin(0.06)
    return c1

def myfunc(tuple_value):
    return tuple_value[1]

def squared_mass_OnlyH(vl) :
    m_H=125.1
    if len(vl)!=2:
        print 'expect 2 vectors exactly for mass calculation for 2H hypothesis, something could be wrong'
        return -1
    else:
        vector_squared_sum=(vl[0].M()-m_H)**2+(vl[1].M()-m_H)**2
        return vector_squared_sum


def squared_mass(vl) :
    m_H=125.1
    m_Z=91.188
    if len(vl)!=3:
        print 'expect 3 vectors exactly for mass calculation, something could be wrong'
        return -1
    else:
        vector_squared_sum=(vl[0].M()-m_H)**2+(vl[1].M()-m_H)**2+(vl[2].M()-m_Z)**2
        vector_squared_sum=min(vector_squared_sum,(vl[0].M()-m_Z)**2+(vl[1].M()-m_H)**2+(vl[2].M()-m_H)**2)
        vector_squared_sum=min(vector_squared_sum,(vl[0].M()-m_H)**2+(vl[1].M()-m_Z)**2+(vl[2].M()-m_H)**2)
        return vector_squared_sum
 

#function which combines 6 jets pairwise into 3 jets, investigate if rechecking jets for better combinations of 5 down to 4 jets
#use BTag information to boost certain jets for H jet combination, to include at least one high Btagged jet
#seems to work less good though than ignoring BTag information, so sideline this function as backup    
def orderedThreeVectorBTagSort(vl_set) :
    result_vector=[]
    index_vector=[]
    case_ind=-1
    m_H=125.1
    m_Z=91.188
    minrange=-1
    vl,BTagList=map(list,zip(*vl_set))

    tuple_value=((),)
    limit_BTag_low=0.60
    if len(vl_set)<6:
        if len(vl_set)>2:
            #print 'expect at least 6 vectors, something needs to be fixed'
            result_vector=[]   
            result_vector.append(vl[0])
            result_vector.append(vl[1])
            result_vector.append(vl[2])   
            index_vector=[0,-1,1,-1,2,-1]
        tuple_value = [(result_vector[0],index_vector[0],index_vector[1]),(result_vector[1],index_vector[2],index_vector[3]),(result_vector[2],index_vector[4],index_vector[5])]
        #print 'before sorting 0 or 5 running',tuple_value[0][0].M(),tuple_value[0][1],tuple_value[0][2],tuple_value[1][0].M(),tuple_value[1][1],tuple_value[1][2],tuple_value[2][0].M(),tuple_value[2][1],tuple_value[2][2]
        tuple_value.sort(key=lambda x:x[0].M(), reverse=True)
        return tuple_value,index_vector
        print 'expect at least 6 vectors, something needs to be fixed'
        return None;
    #either no single high BTag value,or up to 5 high BTag values
    if(vl_set[0][1]<limit_BTag_low or vl_set[4][1]>limit_BTag_low):
        #print 'begin stuff for very high or very low btag',vl_set[0][1],vl_set[4][1]
        perm_list = list(permutations(range(0, 6)))
        case_ind=-1
        for p in perm_list:
            vect1=vl[p[0]]+vl[p[1]]
            vect2=vl[p[2]]+vl[p[3]]
            vect3=vl[p[4]]+vl[p[5]]
                #print 'perm in reco vector ',p,vect1.M(),vect2.M(),vect3.M()
            temp_list=[]   
            temp_list.append(vect1)
            temp_list.append(vect2)
            temp_list.append(vect3)                
            if(minrange<0 or squared_mass(temp_list)<minrange):
                result_vector =  temp_list[:]
                minrange=squared_mass(temp_list)
                index_vector=[]
                index_vector=p
        tuple_value = [(result_vector[0],index_vector[0],index_vector[1]),(result_vector[1],index_vector[2],index_vector[3]),(result_vector[2],index_vector[4],index_vector[5])]
        #print 'before sorting 0 or 5 running',tuple_value[0][0].M(),tuple_value[0][1],tuple_value[0][2],tuple_value[1][0].M(),tuple_value[1][1],tuple_value[1][2],tuple_value[2][0].M(),tuple_value[2][1],tuple_value[2][2]
        tuple_value.sort(key=lambda x:x[0].M(), reverse=True)
        #since here order is not dependent on any BTag info, the H jets are those with the larger H mass
        #print 'after end of 0 or 5 running',tuple_value[0][0].M(),tuple_value[0][1],tuple_value[0][2],tuple_value[1][0].M(),tuple_value[1][1],tuple_value[1][2],tuple_value[2][0].M(),tuple_value[2][1],tuple_value[2][2]
        #HERE HERE HERE HERE HERE study values all of them
        #return tuple_value,case_ind
        if(abs(tuple_value[0][0].M()-m_H)>75 or abs(tuple_value[1][0].M()-m_H)>75):
            if vl_set[0][1]<limit_BTag_low:
                #print 'abnormal high values, case 0',vl_set[0][1],abs(tuple_value[0][0].M()-m_H),abs(tuple_value[1][0].M()-m_H)
                case_ind=-1
            else:
                #print 'abnormal high values, case 1 ',vl_set[0][1],vl_set[1][1],vl_set[2][1],vl_set[3][1],vl_set[4][1],vl_set[5][1],abs(tuple_value[0][0].M()-m_H),abs(tuple_value[1][0].M()-m_H),'single masses',vl_set[0][0].M(),vl_set[1][0].M(),vl_set[2][0].M(),vl_set[3][0].M(),vl_set[4][0].M(),vl_set[5][0].M()
                if vl_set[5][1]<limit_BTag_low:
                    perm_list = list(permutations(range(0, 5)))
                    case_ind=-1
                    for p in perm_list:
                        vect1=vl[p[0]]
                        vect2=vl[p[1]]+vl[p[2]]
                        vect3=vl[p[3]]+vl[p[4]]
                        #print 'perm in reco vector ',p,vect1.M(),vect2.M(),vect3.M()
                        temp_list=[]   
                        temp_list.append(vect1)
                        temp_list.append(vect2)
                        temp_list.append(vect3)                
                        if squared_mass(temp_list)<minrange:
                            minrange=squared_mass(temp_list)
                            result_vector =  temp_list[:]
                            index_vector=[]
                            index_vector=p
                            case_ind=1
                        vect1=vl[p[0]]+vl[p[1]]
                        vect2=vl[p[2]]
                        vect3=vl[p[3]]+vl[p[4]]
                        #print 'perm in reco vector ',p,vect1.M(),vect2.M(),vect3.M()
                        temp_list=[]   
                        temp_list.append(vect1)
                        temp_list.append(vect2)
                        temp_list.append(vect3)                
                        if squared_mass(temp_list)<minrange:
                            minrange=squared_mass(temp_list)
                            result_vector =  temp_list[:]
                            index_vector=[]
                            index_vector=p
                            case_ind=2
                        vect1=vl[p[0]]+vl[p[1]]
                        vect2=vl[p[2]]+vl[p[3]]
                        vect3=vl[p[4]]
                        temp_list=[]   
                        temp_list.append(vect1)
                        temp_list.append(vect2)
                        temp_list.append(vect3)                
                        if squared_mass(temp_list)<minrange:
                            minrange=squared_mass(temp_list)
                            result_vector =  temp_list[:]
                            index_vector=[]
                            index_vector=p
                            case_ind=3
                    if case_ind==1 :
                        tuple_value = [(result_vector[0],index_vector[0],-1),(result_vector[1],index_vector[1],index_vector[2]),(result_vector[2],index_vector[3],index_vector[4])]
                    elif case_ind==2:
                        tuple_value = [(result_vector[0],index_vector[0],index_vector[1]),(result_vector[1],index_vector[2],-1),(result_vector[2],index_vector[3],index_vector[4])]
                    elif case_ind==3:
                        tuple_value = [(result_vector[0],index_vector[0],index_vector[1]),(result_vector[1],index_vector[2],index_vector[3]),(result_vector[2],index_vector[4],-1)]
                    tuple_value.sort(key=lambda x:x[0].M(), reverse=True)
                    #print 'abnormal high values, case 1 but low last value, get to case ',case_ind, index_vector,tuple_value[0][0].M(),tuple_value[1][0].M(),tuple_value[2][0].M()
                #else:
                #    print 'abnormal high values, case 1 but high last value ',vl_set[0][1],vl_set[1][1],vl_set[2][1],vl_set[3][1],vl_set[4][1],vl_set[5][1],abs(tuple_value[0][0].M()-m_H),abs(tuple_value[1][0].M()-m_H)
        return tuple_value,case_ind
    elif(vl_set[0][1]>limit_BTag_low and vl_set[1][1]<limit_BTag_low):
        #only leading jet properly BTagged
        #print 'begin stuff for exactly one high btag',vl_set[0][1],vl_set[1][1]
        perm_list = list(permutations(range(1, 6)))
        for p in perm_list:
            vect1=vl[0]+vl[p[0]]
            vect2=vl[p[1]]+vl[p[2]]
            vect3=vl[p[3]]+vl[p[4]]
            temp_list=[]   
            temp_list.append(vect1)
            mass_sum1=(vect1.M()-m_H)**2+(vect2.M()-m_H)**2+(vect3.M()-m_Z)**2
            mass_sum2=(vect1.M()-m_H)**2+(vect3.M()-m_H)**2+(vect2.M()-m_Z)**2
            if(minrange<0 or min(mass_sum1,mass_sum2)<minrange):
                if(mass_sum1<mass_sum2):
                    temp_list.append(vect2)
                    temp_list.append(vect3)
                    index_vector=[0,p[0],p[1],p[2],p[3],p[4]]
                    minrange=mass_sum1
                else:
                    temp_list.append(vect3)
                    temp_list.append(vect2)
                    index_vector=[0,p[0],p[3],p[4],p[1],p[2]]
                    minrange=mass_sum2
                if minrange!=min(mass_sum1,mass_sum2):
                    print 'wtf for minrange',minrange,min(mass_sum1,mass_sum2)
                result_vector =  temp_list[:]
        tuple_value = [(result_vector[0],index_vector[0],index_vector[1]),(result_vector[1],index_vector[2],index_vector[3]),(result_vector[2],index_vector[4],index_vector[5])]
        #print 'after end of only jet 1 running',tuple_value[0][0].M(),tuple_value[0][1],tuple_value[0][2],tuple_value[1][0].M(),tuple_value[1][1],tuple_value[1][2],tuple_value[2][0].M(),tuple_value[2][1],tuple_value[2][2],vl_set[0][0].M()
        return tuple_value,case_ind
    elif(vl_set[1][1]>limit_BTag_low and vl_set[2][1]<limit_BTag_low):
        #two leading jets properly BTagged
        #print 'begin stuff for exactly two high btags',vl_set[1][1],vl_set[2][1]
        perm_list = list(permutations(range(2, 6)))
        for p in perm_list:
            temp_list=[]   
            vect1=vl[0]+vl[p[0]]
            vect2=vl[1]+vl[p[1]]
            vect3=vl[p[2]]+vl[p[3]]
            mass_sum=(vect1.M()-m_H)**2+(vect2.M()-m_H)**2+(vect3.M()-m_Z)**2
            #print 'indices checked 0 here',[0,p[0],1,p[1],p[2],p[3]],vect1.M(),vect2.M(),vect3.M()
            if(minrange<0 or mass_sum<minrange):
                temp_list.append(vect1)
                temp_list.append(vect2)
                temp_list.append(vect3)
                index_vector=[0,p[0],1,p[1],p[2],p[3]]
                minrange=mass_sum
                result_vector =  temp_list[:]
        if(result_vector[0].M()>result_vector[1].M()):
            tuple_value = [(result_vector[0],index_vector[0],index_vector[1]),(result_vector[1],index_vector[2],index_vector[3]),(result_vector[2],index_vector[4],index_vector[5])]
        else:
            tuple_value = [(result_vector[1],index_vector[2],index_vector[3]),(result_vector[0],index_vector[0],index_vector[1]),(result_vector[2],index_vector[4],index_vector[5])]
        #HERE HERE HERE HERE HERE study values all of them
        #return tuple_value,case_ind
        if(abs(tuple_value[0][0].M()-m_H)>75 or abs(tuple_value[1][0].M()-m_H)>75) or abs(tuple_value[2][0].M()-m_Z)>75:
            #print 'abnormal high values for higgs or Z, 2jets',vl_set[0][0].M(),vl_set[0][0].E(),vl_set[1][0].M(),vl_set[1][0].E(),vl_set[2][0].M(),vl_set[1][0].E(),vl_set[3][0].M(),vl_set[3][0].E(),vl_set[4][0].M(),vl_set[5][0].E(),vl_set[5][0].M(),vl_set[5][0].E(),abs(tuple_value[0][0].M()-m_H),abs(tuple_value[1][0].M()-m_H),abs(tuple_value[2][0].M()-m_Z),(vl_set[0][0]+vl_set[1][0]+vl_set[2][0]+vl_set[3][0]+vl_set[4][0]+vl_set[5][0]).M()
            #solution -> cut out one of the jets from the combinations
            perm_list = list(permutations(range(2, 6)))
            case_ind=-1
            for p in perm_list:
                temp_list=[]   
                vect1=vl[0]+vl[p[0]]
                vect2=vl[1]+vl[p[1]]
                vect3=vl[p[2]]
                mass_sum=(vect1.M()-m_H)**2+(vect2.M()-m_H)**2+(vect3.M()-m_Z)**2
                #print 'indices checked 2 30 here get',[0,p[0],1,p[1],p[2]],vect1.M(),vect2.M(),vect3.M(), mass_sum,minrange
                if mass_sum<minrange:
                    temp_list.append(vect1)
                    temp_list.append(vect2)
                    temp_list.append(vect3)
                    case_ind=3
                    index_vector=[0,p[0],1,p[1],p[2],-1]
                    result_vector =  temp_list[:]
                    minrange=mass_sum
                    #print 'indices checked 2 30 here',[0,p[0],1,p[1],p[2]],vect1.M(),vect2.M(),vect3.M(), 'new min',minrange
                vect1=vl[0]+vl[p[0]]
                vect2=vl[1]+vl[p[1]]
                vect3=vl[p[3]]
                mass_sum=(vect1.M()-m_H)**2+(vect2.M()-m_H)**2+(vect3.M()-m_Z)**2
                #print 'indices checked 2 31 here get',[0,p[0],1,p[1],p[3]],vect1.M(),vect2.M(),vect3.M(), mass_sum,minrange
                if mass_sum<minrange:
                    temp_list.append(vect1)
                    temp_list.append(vect2)
                    temp_list.append(vect3)
                    case_ind=3
                    index_vector=[0,p[0],1,p[1],p[3],-1]
                    result_vector =  temp_list[:]
                    minrange=mass_sum
                    #print 'indices checked 2 31 here',[0,p[0],1,p[1],p[3]],vect1.M(),vect2.M(),vect3.M(), 'new min',minrange
                vect1=vl[0]+vl[p[0]]
                vect2=vl[1]
                vect3=vl[p[2]]+vl[p[3]]
                mass_sum=(vect1.M()-m_H)**2+(vect2.M()-m_H)**2+(vect3.M()-m_Z)**2
                #print 'indices checked 2 20 here get',[0,p[0],1,p[2],p[3]],vect1.M(),vect2.M(),vect3.M(), mass_sum,minrange
                if mass_sum<minrange:
                    temp_list.append(vect1)
                    temp_list.append(vect2)
                    temp_list.append(vect3)
                    case_ind=2
                    index_vector=[0,p[0],1,-1,p[2],p[3]]
                    result_vector =  temp_list[:]
                    minrange=mass_sum
                    #print 'indices checked 2 20 here',[0,p[0],1,p[2],p[3]],vect1.M(),vect2.M(),vect3.M(), 'new min',minrange
                vect1=vl[0]
                vect2=vl[1]+vl[p[1]]
                vect3=vl[p[2]]+vl[p[3]]
                mass_sum=(vect1.M()-m_H)**2+(vect2.M()-m_H)**2+(vect3.M()-m_Z)**2
                #print 'indices checked 2 10 here get',[0,1,p[1],p[2],p[3]],vect1.M(),vect2.M(),vect3.M(), mass_sum,minrange
                if mass_sum<minrange:
                    temp_list.append(vect1)
                    temp_list.append(vect2)
                    temp_list.append(vect3)
                    case_ind=1
                    index_vector=[0,-1,1,p[1],p[2],p[3]]
                    result_vector =  temp_list[:]
                    minrange=mass_sum
                    #print 'indices checked 2 10 here',[0,1,p[1],p[2],p[3]],vect1.M(),vect2.M(),vect3.M(), 'new min',minrange
            #index vectors properly filled with -1 already at that point
            if(result_vector[0].M()>result_vector[1].M()):
                tuple_value = [(result_vector[0],index_vector[0],index_vector[1]),(result_vector[1],index_vector[2],index_vector[3]),(result_vector[2],index_vector[4],index_vector[5])]
            else:
                tuple_value = [(result_vector[1],index_vector[2],index_vector[3]),(result_vector[0],index_vector[0],index_vector[1]),(result_vector[2],index_vector[4],index_vector[5])]
            #print 'after end of exactly jet 2 running with reschedule',tuple_value[0][0].M(),tuple_value[1][0].M(),tuple_value[2][0].M(),(tuple_value[0][0]+tuple_value[1][0]+tuple_value[2][0]).M(),"orig",(vl_set[0][0]+vl_set[1][0]+vl_set[2][0]+vl_set[3][0]+vl_set[4][0]+vl_set[5][0]).M()
        return tuple_value,case_ind
    elif vl_set[2][1]>limit_BTag_low :
        #three leading jets properly BTagged, but also as maximum 4 b-jets, leaves possibility of one 0 b-tagged combination
        #print 'begin stuff for three or four btags',vl_set[0][1],vl_set[1][1],vl_set[2][1],vl_set[3][1]
        if vl_set[4][1]>limit_BTag_low :
            print 'WRONG WRONG WRONG, should have been taken care of previously'
        case_ind=-1
        perm_list = list(permutations(range(0, 6)))
        for p in perm_list:
            temp_list=[]   
            vect1=vl[p[0]]+vl[p[1]]
            vect2=vl[p[2]]+vl[p[3]]
            vect3=vl[p[4]]+vl[p[5]]
            case_ind=-1
            #check if any combination has non b-tag values
            if(vl_set[p[0]][1]<limit_BTag_low and vl_set[p[1]][1]<limit_BTag_low):
                mass_sum=(vect2.M()-m_H)**2+(vect3.M()-m_H)**2+(vect1.M()-m_Z)**2 
                if(minrange<0 or mass_sum<minrange):
                    min_range=mass_sum
                    if vect2.M()>vect3.M():
                        temp_list.append(vect2)
                        temp_list.append(vect3)
                        temp_list.append(vect1)
                        index_vector=[p[2],p[3],p[4],p[5],p[0],p[1]]
                        result_vector =  temp_list[:]
                    else: 
                        temp_list.append(vect3)
                        temp_list.append(vect2)
                        temp_list.append(vect1)
                        index_vector=[p[4],p[5],p[2],p[3],p[0],p[1]]
                        result_vector =  temp_list[:]
            elif (vl_set[p[2]][1]<limit_BTag_low and vl_set[p[3]][1]<limit_BTag_low):
                mass_sum=(vect1.M()-m_H)**2+(vect3.M()-m_H)**2+(vect2.M()-m_Z)**2 
                if(minrange<0 or mass_sum<minrange):
                    min_range=mass_sum
                    if vect1.M()>vect3.M():
                        temp_list.append(vect1)
                        temp_list.append(vect3)
                        temp_list.append(vect2)
                        index_vector=[p[0],p[1],p[4],p[5],p[2],p[3]]
                        result_vector =  temp_list[:]
                    else: 
                        temp_list.append(vect3)
                        temp_list.append(vect1)
                        temp_list.append(vect2)
                        index_vector=[p[4],p[5],p[0],p[1],p[2],p[3]]
                        result_vector =  temp_list[:]
            elif (vl_set[p[4]][1]<limit_BTag_low and vl_set[p[5]][1]<limit_BTag_low):
                mass_sum=(vect1.M()-m_H)**2+(vect2.M()-m_H)**2+(vect3.M()-m_Z)**2 
                if(minrange<0 or mass_sum<minrange):
                    min_range=mass_sum
                    if vect1.M()>vect2.M():
                        temp_list.append(vect1)
                        temp_list.append(vect2)
                        temp_list.append(vect3)
                        index_vector=[p[0],p[1],p[2],p[3],p[4],p[5]]
                        result_vector =  temp_list[:]
                    else: 
                        temp_list.append(vect2)
                        temp_list.append(vect1)
                        temp_list.append(vect3)
                        index_vector=[p[1],p[2],p[0],p[1],p[4],p[5]]
                        result_vector =  temp_list[:]
            else: 
                #here all jets have at least a b-jet
                temp_list=[]   
                temp_list.append(vect1)
                temp_list.append(vect2)
                temp_list.append(vect3)                
                if(minrange<0 or squared_mass(temp_list)<minrange):
                    minrange=squared_mass(temp_list)
                    result_vector =  temp_list[:]
                    result_vector.sort(key=lambda x: x.M(), reverse=True)
                    #print'resultvec', result_vector[0].M(),result_vector[1].M(),result_vector[2].M()
                    if vect1.M()>vect2.M() and vect1.M()>vect3.M():
                        if(vect2.M()>vect3.M()):
                            index_vector=[p[0],p[1],p[2],p[3],p[4],p[5]]
                        else:
                            index_vector=[p[0],p[1],p[4],p[5],p[2],p[3]]
                    elif vect2.M()>vect1.M() and vect2.M()>vect3.M():
                        if(vect1.M()>vect3.M()):
                            index_vector=[p[2],p[3],p[0],p[1],p[4],p[5]]
                        else:
                            index_vector=[p[2],p[3],p[4],p[5],p[0],p[1]]
                    elif vect3.M()>vect1.M() and vect3.M()>vect2.M():
                       if(vect1.M()>vect2.M()):
                            index_vector=[p[4],p[5],p[0],p[1],p[2],p[3]]
                       else:
                            index_vector=[p[4],p[5],p[2],p[3],p[0],p[1]]
                    else:
                        print 'what went wrong with the index order',vect1.M(),vect2.M(),vect3.M()
        tuple_value = [(result_vector[0],index_vector[0],index_vector[1]),(result_vector[1],index_vector[2],index_vector[3]),(result_vector[2],index_vector[4],index_vector[5])]
        return tuple_value,case_ind
    if len(result_vector)<3 or len(index_vector)<5:
        print 'what the hell is wrong'

#works better as default, combine n jets with n>=4 into 3 jets, assume default is 6 jets
#ignore BTag information but keep track of indices which end up in the jet mass combinations
#order jets according to the masses of the combined jets
def orderedThreeVector(vl) :
    result_vector=[]
    index_vector=[]
    m_H=125.1
    m_Z=91.188
    if len(vl)<3:
        print 'expect at least 3 vectors, something needs to be fixed'
        return None;
    if len(vl)==3:
        result_vector = vl
        index_vector= [0,-1,1,-1,2,-1]
    else:
        if len(vl)==4:
            minrange=-1
            for vect in (range(len(vl))):
                #print 'in index',vect,'of',len(vl)
                for vect2 in range(vect+1,len(vl),1): 
                    temp_list=vl[:]
                    #print 'vect,vect2 at start',vect,vect2,len(vl),len(temp_list)
                    temp_list[vect]=temp_list[vect]+temp_list[vect2]
                    temp_list.pop(vect2)      
                    index_vector.append(-1)
                    #print 'vect,vect2',vect,vect2,len(temp_list)
                #for vect in len(range(vl-1)):
                    index_vector=[vect,vect2]
                    if(minrange<0):
                        result_vector =  temp_list
                        #print 'get in here once',vect,vect2,squared_mass(temp_list),minrange
                        minrange=squared_mass(temp_list)
                        for vect_temp in (range(len(vl))):
                            if(vect_temp!=vect and vect_temp!=vect2):
                                index_vector.append(vect_temp)
                        #print 'get in here once 2',vect,vect2,squared_mass(temp_list),minrange
                    elif squared_mass(temp_list)<minrange:
                        #print 'get in here maybe',vect,vect2,squared_mass(temp_list),minrange
                        minrange=squared_mass(temp_list)
                        result_vector =  temp_list
                        for vect_temp in (range(len(vl))):
                            if(vect_temp!=vect and vect_temp!=vect2):
                                index_vector.append(vect_temp)
            index_vector= [0,-1,1,-1,2,-1]
        elif len(vl)==5:
            minrange=-1
            #print 'begin stuff here'
            for vect in (range(len(vl))):
                #print 'in index',vect,'of',len(vl)
                for vect2 in range(vect+1,len(vl),1): 
                    #print 'zero check',vect,vect+1,vect2
                    #print 'vect,vect2 at start',vect,vect2,len(vl)
                    for vect3 in range(vect+1,len(vl),1): 
                        #print 'first draft',vect,vect2,vect3
                        if vect3!=vect2 :
                            #print 'first check',vect,vect2,vect3
                            for vect4 in range(vect3+1,len(vl)):
                                #print 'second draft',vect,vect2,vect3,vect4
                                if vect4!=vect2:
                                    #print 'vect,vect2 at start',vect,vect2,len(vl),len(temp_list)
                                    temp_list=vl[:]
                                    #print 'second check',vect,vect2,vect3,vect4,len(vl),len(temp_list)
                                    temp_list[vect]=temp_list[vect]+temp_list[vect2]
                                    temp_list[vect3]=temp_list[vect3]+temp_list[vect4]
                                    if(vect2<vect4):
                                        #print 'pop1',vect,vect2,vect3,vect4,len(temp_list)
                                        temp_list.pop(vect4)      
                                        #print 'pop1 vect4 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                        temp_list.pop(vect2) 
                                        #print 'pop1 vect2 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                    else:
                                        #print 'pop2',vect,vect2,vect3,vect4,len(temp_list)
                                        temp_list.pop(vect2)      
                                        #print 'pop2 vect2 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                        temp_list.pop(vect4)  
                                        #print 'pop2 vect4 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                    if(minrange<0):
                                        result_vector =  temp_list
                                        #print 'get in here once',vect,vect2,squared_mass(temp_list),minrange
                                        minrange=squared_mass(temp_list)
                                        index_vector=[]
                                        index_vector=[vect,vect2,vect3,vect4]
                                        for vect_temp in (range(len(vl))):
                                            if(vect_temp!=vect and vect_temp!=vect2 and vect_temp!=vect3 and vect_temp!=vect4):
                                                index_vector.append(vect_temp)
                                        #print 'get in here once 2',vect,vect2,squared_mass(temp_list),minrange
                                    elif squared_mass(temp_list)<minrange:
                                        #print 'get in here maybe',vect,vect2,squared_mass(temp_list),minrange
                                        for vect_temp in (range(len(vl))):
                                            if(vect_temp!=vect and vect_temp!=vect2 and vect_temp!=vect3 and vect_temp!=vect4):
                                                index_vector.append(vect_temp)
                                        minrange=squared_mass(temp_list)
                                        result_vector =  temp_list
            index_vector.append(-1)
        elif len(vl)==6:
            minrange=-1
            #print 'begin stuff here'
            perm_list = list(permutations(range(0, 6)))
            for p in perm_list:
                vect1=vl[p[0]]+vl[p[1]]
                vect2=vl[p[2]]+vl[p[3]]
                vect3=vl[p[4]]+vl[p[5]]
                #print 'perm in reco vector ',p,vect1.M(),vect2.M(),vect3.M()
                temp_list=[]   
                temp_list.append(vect1)
                temp_list.append(vect2)
                temp_list.append(vect3)                
                if(minrange<0):
                    result_vector =  temp_list
                    minrange=squared_mass(temp_list)
                    index_vector=[]
                    index_vector=p
                elif squared_mass(temp_list)<minrange:
                    minrange=squared_mass(temp_list)
                    result_vector =  temp_list
                    index_vector=[]
                    index_vector=p

    tuple_value = [(result_vector[0],index_vector[0],index_vector[1]),(result_vector[1],index_vector[2],index_vector[3]),(result_vector[2],index_vector[4],index_vector[5])]
    tuple_value.sort(key=lambda x:x[0].M(), reverse=True)
    return tuple_value
 
def CLICdpStyle(): 
    gROOT.SetStyle("Plain") 
    gStyle.SetCanvasColor(root.kWhite)
    gStyle.SetFrameFillColor(root.kWhite)
    gStyle.SetStatColor(root.kWhite)
    gStyle.SetPadColor(root.kWhite)
    gStyle.SetFillColor(10)
    gStyle.SetTitleFillColor(root.kWhite)
  
    gStyle.SetPaperSize(20, 26) 
    
    gStyle.SetDrawBorder(0)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetFrameBorderMode(0)
    gStyle.SetLegendBorderSize(0)
  
    gStyle.SetTextSize(0.05)
    gStyle.SetTitleSize(0.06,"xyz")
    gStyle.SetLabelSize(0.06,"xyz")
    
    gStyle.SetLabelOffset(0.015,"xyz")
    gStyle.SetTitleOffset(1.2,"yz") 
    gStyle.SetTitleOffset(1.17,"x")
 
    font = 42 
    gStyle.SetTitleFont(font)
    gStyle.SetTitleFontSize(0.06)
    gStyle.SetStatFont(font)
    gStyle.SetStatFontSize(0.07)
    gStyle.SetTextFont(font)
    gStyle.SetLabelFont(font,"xyz")
    gStyle.SetTitleFont(font,"xyz")
    gStyle.SetTitleBorderSize(0)
    gStyle.SetStatBorderSize(1)
    gStyle.SetMarkerStyle(1)
    gStyle.SetLineWidth(2)  
    gStyle.SetMarkerSize(1.2)
    gStyle.SetPalette(1) 

    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0) 
    gStyle.SetOptFit(0) 
    gStyle.SetEndErrorSize(5)   

    gStyle.SetHistLineWidth(2)
    gStyle.SetFrameLineWidth(2)
    gStyle.SetFuncWidth(2)
    gStyle.SetHistLineColor(root.kBlack)
    gStyle.SetFuncColor(root.kBlack)
    gStyle.SetLabelColor(root.kBlack,"xyz")


    gStyle.SetPadBottomMargin(0.18)
    gStyle.SetPadTopMargin(0.11)
    gStyle.SetPadRightMargin(0.08)
    gStyle.SetPadLeftMargin(0.17)
   
    
    gStyle.SetNdivisions(506, "xy")
   
    gStyle.SetPadGridX(0)
    gStyle.SetPadGridY(0)
   
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    gStyle.SetCanvasDefW(800)
    gStyle.SetCanvasDefH(700)
    
    gROOT.ForceStyle()

def fill_HHZ_histograms(file,xsec,mytree,hist_vec_reco_1D_parton,lumi,ishhzfile):
#, h_hist_parton, h_hist_vec_gen, h_hist_vec_reco_1D_parton, usePartonInfo, x_sec, fill_partonInfo, fill_genInfo):
#    print "size of histogram vectors ",len(h_hist_parton),len(h_hist_vec_gen),len(h_hist_vec_reco_1D_parton), " booleans ",usePartonInfo,fill_partonInfo,fill_genInfo, " xsec ",xsec 
    
    use_EMissNeutrinoProjection=True
    #here use total 4 vector, or 4 vector sum of jet 1 and 2 (see flag below)- isolated photon four vector plus correction with EMiss on both jet axes
    #mass cuts are then also done after projecting the EMiss
    use_sqrtJets=True #in this case use j1+j2-isolated photons and with upper flag still decide if EMiss projection on jets is performed
    use_MHMiss_over_PFOMiss = False #if recoil jet missing energy or PFO missing energy is used in the missing energy projection

    fCut_mass_1_min=0.
    fCut_mass1_center=126.
    fCut_mass2_center=92.5
    fCut_mass1_radius=35.
    fCut_mass2_radius=35.

    fCut_thetaWindow=70.
    fCut_thetaRef=90.
    fCut_delta_theta = 100.
    
    t_var_weight  = array('f',[0])
    t_var_sqrtS  = array('f',[0])
    t_var_sqrtS_jets  = array('f',[0])
    t_var_sqrtS_gen  = array('f',[0])
    t_var_sqrtS_genjets  = array('f',[0])
    t_var_sqrtS_parton  = array('f',[0])
    t_var_MET = array('f',[0])
    t_var_MET_gen = array('f',[0])
    t_var_MHT = array('f',[0])
    t_var_MHT_gen = array('f',[0])

    t_var_comb_jet1_mass = array('f',[0])
    t_var_comb_jet1_theta = array('f',[0])
    t_var_comb_jet1_phi = array('f',[0])
    t_var_comb_jet1_E = array('f',[0])
    t_var_comb_jet1_BTagMax = array('f',[0])
    #angles between original jets combined in this one
    t_var_comb_jet1_dalpha = array('f',[0])
    t_var_comb_jet1_dphi = array('f',[0])
    t_var_comb_jet1_dtheta = array('f',[0])
    #ratio of more energetic input jet for combination to total sum
    t_var_comb_jet1_E1_over_Etot = array('f',[0])
    #t_var_comb_jet1_jetChargeMax = array('f',[0])
    #t_var_comb_jet1_jetChargeMin = array('f',[0])
    #t_var_comb_jet1_E1_jetCharge = array('f',[0])

    t_var_comb_jet2_mass = array('f',[0])
    t_var_comb_jet2_theta = array('f',[0])
    t_var_comb_jet2_phi = array('f',[0])
    t_var_comb_jet2_E = array('f',[0])
    t_var_comb_jet2_BTagMax = array('f',[0])
    t_var_comb_jet2_dalpha = array('f',[0])
    t_var_comb_jet2_dphi = array('f',[0])
    t_var_comb_jet2_dtheta = array('f',[0])
    t_var_comb_jet2_E1_over_Etot = array('f',[0])
    #t_var_comb_jet2_jetChargeMax = array('f',[0])
    #t_var_comb_jet2_jetChargeMin = array('f',[0])
    #t_var_comb_jet2_E1_jetCharge = array('f',[0])

    t_var_comb_jet3_mass = array('f',[0])
    t_var_comb_jet3_theta = array('f',[0])
    t_var_comb_jet3_phi = array('f',[0])
    t_var_comb_jet3_E = array('f',[0])
    t_var_comb_jet3_BTagMax = array('f',[0])
    t_var_comb_jet3_dalpha = array('f',[0])
    t_var_comb_jet3_dphi = array('f',[0])
    t_var_comb_jet3_dtheta = array('f',[0])
    t_var_comb_jet3_E1_over_Etot = array('f',[0])
    #t_var_comb_jet3_jetChargeMax = array('f',[0])
    #t_var_comb_jet3_jetChargeMin = array('f',[0])
    #t_var_comb_jet3_E1_jetCharge = array('f',[0])

    t_var_comb_genjet1_mass = array('f',[0])
    t_var_comb_genjet1_theta = array('f',[0])
    t_var_comb_genjet1_phi = array('f',[0])
    t_var_comb_genjet1_E = array('f',[0])
    #angles between original jets combined in this one
    t_var_comb_genjet1_dalpha = array('f',[0])
    t_var_comb_genjet1_dphi = array('f',[0])
    t_var_comb_genjet1_dtheta = array('f',[0])
    #ratio of more energetic input jet for combination to total sum
    t_var_comb_genjet1_E1_over_Etot = array('f',[0])

    t_var_comb_genjet2_mass = array('f',[0])
    t_var_comb_genjet2_theta = array('f',[0])
    t_var_comb_genjet2_phi = array('f',[0])
    t_var_comb_genjet2_E = array('f',[0])
    #angles between original jets combined in this one
    t_var_comb_genjet2_dalpha = array('f',[0])
    t_var_comb_genjet2_dphi = array('f',[0])
    t_var_comb_genjet2_dtheta = array('f',[0])
    #ratio of more energetic input jet for combination to total sum
    t_var_comb_genjet2_E1_over_Etot = array('f',[0])

    t_var_comb_genjet3_mass = array('f',[0])
    t_var_comb_genjet3_theta = array('f',[0])
    t_var_comb_genjet3_phi = array('f',[0])
    t_var_comb_genjet3_E = array('f',[0])
    #angles between original jets combined in this one
    t_var_comb_genjet3_dalpha = array('f',[0])
    t_var_comb_genjet3_dphi = array('f',[0])
    t_var_comb_genjet3_dtheta = array('f',[0])
    #ratio of more energetic input jet for combination to total sum
    t_var_comb_genjet3_E1_over_Etot = array('f',[0])

    #order jets here by energy
    t_var_jet1_E = array('f',[0])
    t_var_jet1_theta = array('f',[0])
    t_var_jet1_BTag = array('f',[0])
    #t_var_jet1_jetCharge = array('f',[0])

    t_var_jet2_E = array('f',[0])
    t_var_jet2_theta = array('f',[0])
    t_var_jet2_BTag = array('f',[0])
    #t_var_jet2_jetCharge = array('f',[0])

    t_var_jet3_E = array('f',[0])
    t_var_jet3_theta = array('f',[0])
    t_var_jet3_BTag = array('f',[0])
    #t_var_jet3_jetCharge = array('f',[0])

    t_var_jet4_E = array('f',[0])
    t_var_jet4_theta = array('f',[0])
    t_var_jet4_BTag = array('f',[0])
    #t_var_jet4_jetCharge = array('f',[0])

    t_var_jet5_E = array('f',[0])
    t_var_jet5_theta = array('f',[0])
    t_var_jet5_BTag = array('f',[0])
    #t_var_jet5_jetCharge = array('f',[0])

    t_var_jet6_E = array('f',[0])
    t_var_jet6_theta = array('f',[0])
    t_var_jet6_BTag = array('f',[0])
    #t_var_jet6_jetCharge = array('f',[0])

    t_var_y12 = array('f',[0])
    t_var_y23 = array('f',[0])
    t_var_y34 = array('f',[0])
    t_var_y45 = array('f',[0])
    t_var_y56 = array('f',[0])

    t_var_gen_y12 = array('f',[0])
    t_var_gen_y23 = array('f',[0])
    t_var_gen_y34 = array('f',[0])
    t_var_gen_y45 = array('f',[0])
    t_var_gen_y56 = array('f',[0])

    #order by BTag --> largest BTag is called BTag1
    t_var_BTag1 = array('f',[0])
    t_var_BTag2 = array('f',[0])
    t_var_BTag3 = array('f',[0])
    t_var_BTag4 = array('f',[0])

    t_var_BTag_sum_max2 = array('f',[0])
    t_var_BTag_sum_max3 = array('f',[0])
    t_var_BTag_sum_max4 = array('f',[0])
    t_var_BTag_sum_all = array('f',[0])

    mytree.Branch('weight', t_var_weight , 'weight/F')
    mytree.Branch('sqrtS', t_var_sqrtS , 'sqrtS/F')
    mytree.Branch('sqrtS_jets', t_var_sqrtS_jets , 'sqrtS_jets/F')
    mytree.Branch('sqrtS_gen', t_var_sqrtS_gen , 'sqrtS_gen/F')
    mytree.Branch('sqrtS_genjets', t_var_sqrtS_genjets , 'sqrtS_genjets/F')
    mytree.Branch('sqrtS_parton', t_var_sqrtS_parton , 'sqrtS_parton/F')

    mytree.Branch('MET',t_var_MET,'MET/F')
    mytree.Branch('MET_gen',t_var_MET_gen,'MET_gen/F')
    mytree.Branch('MHT',t_var_MHT,'MHT/F')
    mytree.Branch('MHT_gen',t_var_MHT_gen,'MHT_gen/F')

    mytree.Branch('comb_jet1_mass',t_var_comb_jet1_mass,'comb_jet1_mass/F')
    mytree.Branch('comb_jet1_theta',t_var_comb_jet1_theta,'comb_jet1_theta/F')
    mytree.Branch('comb_jet1_phi',t_var_comb_jet1_phi,'comb_jet1_phi/F')
    mytree.Branch('comb_jet1_E',t_var_comb_jet1_E,'comb_jet1_E/F')
    mytree.Branch('comb_jet1_BTagMax',t_var_comb_jet1_BTagMax,'comb_jet1_BTagMax/F')
    mytree.Branch('comb_jet1_dalpha',t_var_comb_jet1_dalpha,'comb_jet1_dalpha/F')
    mytree.Branch('comb_jet1_dphi',t_var_comb_jet1_dphi,'comb_jet1_dphi/F')
    mytree.Branch('comb_jet1_dtheta',t_var_comb_jet1_dtheta,'comb_jet1_dtheta/F')
    mytree.Branch('comb_jet1_E1_over_Etot',t_var_comb_jet1_E1_over_Etot,'comb_jet1_E1_over_Etot/F')
    #mytree.Branch('comb_jet1_jetChargeMax',t_var_comb_jet1_jetChargeMax,'comb_jet1_jetChargeMax/F')
    #mytree.Branch('comb_jet1_jetChargeMin',t_var_comb_jet1_jetChargeMin,'comb_jet1_jetChargeMin/F')
    #mytree.Branch('comb_jet1_E1_jetCharge',t_var_comb_jet1_E1_jetCharge,'comb_jet1_E1_jetCharge/F')

    mytree.Branch('comb_jet2_mass',t_var_comb_jet2_mass,'comb_jet2_mass/F')
    mytree.Branch('comb_jet2_theta',t_var_comb_jet2_theta,'comb_jet2_theta/F')
    mytree.Branch('comb_jet2_phi',t_var_comb_jet2_phi,'comb_jet2_phi/F')
    mytree.Branch('comb_jet2_E',t_var_comb_jet2_E,'comb_jet2_E/F')
    mytree.Branch('comb_jet2_BTagMax',t_var_comb_jet2_BTagMax,'comb_jet2_BTagMax/F')
    mytree.Branch('comb_jet2_dalpha',t_var_comb_jet2_dalpha,'comb_jet2_dalpha/F')
    mytree.Branch('comb_jet2_dphi',t_var_comb_jet2_dphi,'comb_jet2_dphi/F')
    mytree.Branch('comb_jet2_dtheta',t_var_comb_jet2_dtheta,'comb_jet2_dtheta/F')
    mytree.Branch('comb_jet2_E1_over_Etot',t_var_comb_jet2_E1_over_Etot,'comb_jet2_E1_over_Etot/F')
    #mytree.Branch('comb_jet2_jetChargeMax',t_var_comb_jet2_jetChargeMax,'comb_jet2_jetChargeMax/F')
    #mytree.Branch('comb_jet2_jetChargeMin',t_var_comb_jet2_jetChargeMin,'comb_jet2_jetChargeMin/F')
    #mytree.Branch('comb_jet2_E1_jetCharge',t_var_comb_jet2_E1_jetCharge,'comb_jet2_E1_jetCharge/F')

    mytree.Branch('comb_jet3_mass',t_var_comb_jet3_mass,'comb_jet3_mass/F')
    mytree.Branch('comb_jet3_theta',t_var_comb_jet3_theta,'comb_jet3_theta/F')
    mytree.Branch('comb_jet3_phi',t_var_comb_jet3_phi,'comb_jet3_phi/F')
    mytree.Branch('comb_jet3_E',t_var_comb_jet3_E,'comb_jet3_E/F')
    mytree.Branch('comb_jet3_BTagMax',t_var_comb_jet3_BTagMax,'comb_jet3_BTagMax/F')
    mytree.Branch('comb_jet3_dalpha',t_var_comb_jet3_dalpha,'comb_jet3_dalpha/F')
    mytree.Branch('comb_jet3_dphi',t_var_comb_jet3_dphi,'comb_jet3_dphi/F')
    mytree.Branch('comb_jet3_dtheta',t_var_comb_jet3_dtheta,'comb_jet3_dtheta/F')
    mytree.Branch('comb_jet3_E1_over_Etot',t_var_comb_jet3_E1_over_Etot,'comb_jet3_E1_over_Etot/F')
    #mytree.Branch('comb_jet3_jetChargeMax',t_var_comb_jet3_jetChargeMax,'comb_jet3_jetChargeMax/F')
    #mytree.Branch('comb_jet3_jetChargeMin',t_var_comb_jet3_jetChargeMin,'comb_jet3_jetChargeMin/F')
    #mytree.Branch('comb_jet3_E1_jetCharge',t_var_comb_jet3_E1_jetCharge,'comb_jet3_E1_jetCharge/F')

    mytree.Branch('comb_genjet1_mass',t_var_comb_genjet1_mass,'comb_genjet1_mass/F')
    mytree.Branch('comb_genjet1_theta',t_var_comb_genjet1_theta,'comb_genjet1_theta/F')
    mytree.Branch('comb_genjet1_phi',t_var_comb_genjet1_phi,'comb_genjet1_phi/F')
    mytree.Branch('comb_genjet1_E',t_var_comb_genjet1_E,'comb_genjet1_E/F')
    mytree.Branch('comb_genjet1_dalpha',t_var_comb_genjet1_dalpha,'comb_genjet1_dalpha/F')
    mytree.Branch('comb_genjet1_dphi',t_var_comb_genjet1_dphi,'comb_genjet1_dphi/F')
    mytree.Branch('comb_genjet1_dtheta',t_var_comb_genjet1_dtheta,'comb_genjet1_dtheta/F')
    mytree.Branch('comb_genjet1_E1_over_Etot',t_var_comb_genjet1_E1_over_Etot,'comb_genjet1_E1_over_Etot/F')

    mytree.Branch('comb_genjet2_mass',t_var_comb_genjet2_mass,'comb_genjet2_mass/F')
    mytree.Branch('comb_genjet2_theta',t_var_comb_genjet2_theta,'comb_genjet2_theta/F')
    mytree.Branch('comb_genjet2_phi',t_var_comb_genjet2_phi,'comb_genjet2_phi/F')
    mytree.Branch('comb_genjet2_E',t_var_comb_genjet2_E,'comb_genjet2_E/F')
    mytree.Branch('comb_genjet2_dalpha',t_var_comb_genjet2_dalpha,'comb_genjet2_dalpha/F')
    mytree.Branch('comb_genjet2_dphi',t_var_comb_genjet2_dphi,'comb_genjet2_dphi/F')
    mytree.Branch('comb_genjet2_dtheta',t_var_comb_genjet2_dtheta,'comb_genjet2_dtheta/F')
    mytree.Branch('comb_genjet2_E1_over_Etot',t_var_comb_genjet2_E1_over_Etot,'comb_genjet2_E1_over_Etot/F')

    mytree.Branch('comb_genjet3_mass',t_var_comb_genjet3_mass,'comb_genjet3_mass/F')
    mytree.Branch('comb_genjet3_theta',t_var_comb_genjet3_theta,'comb_genjet3_theta/F')
    mytree.Branch('comb_genjet3_phi',t_var_comb_genjet3_phi,'comb_genjet3_phi/F')
    mytree.Branch('comb_genjet3_E',t_var_comb_genjet3_E,'comb_genjet3_E/F')
    mytree.Branch('comb_genjet3_dalpha',t_var_comb_genjet3_dalpha,'comb_genjet3_dalpha/F')
    mytree.Branch('comb_genjet3_dphi',t_var_comb_genjet3_dphi,'comb_genjet3_dphi/F')
    mytree.Branch('comb_genjet3_dtheta',t_var_comb_genjet3_dtheta,'comb_genjet3_dtheta/F')
    mytree.Branch('comb_genjet3_E1_over_Etot',t_var_comb_genjet3_E1_over_Etot,'comb_genjet3_E1_over_Etot/F')

    mytree.Branch('jet1_E',t_var_jet1_E,'jet1_E/F')
    mytree.Branch('jet1_theta',t_var_jet1_theta,'jet1_theta/F')
    mytree.Branch('jet1_BTag',t_var_jet1_BTag,'jet1_BTag/F')
    #mytree.Branch('jet1_jetCharge',t_var_jet1_jetCharge,'jet1_jetCharge/F')

    mytree.Branch('jet2_E',t_var_jet2_E,'jet2_E/F')
    mytree.Branch('jet2_theta',t_var_jet2_theta,'jet2_theta/F')
    mytree.Branch('jet2_BTag',t_var_jet2_BTag,'jet2_BTag/F')
    #mytree.Branch('jet2_jetCharge',t_var_jet2_jetCharge,'jet2_jetCharge/F')

    mytree.Branch('jet3_E',t_var_jet3_E,'jet3_E/F')
    mytree.Branch('jet3_theta',t_var_jet3_theta,'jet3_theta/F')
    mytree.Branch('jet3_BTag',t_var_jet3_BTag,'jet3_BTag/F')
    #mytree.Branch('jet3_jetCharge',t_var_jet3_jetCharge,'jet3_jetCharge/F')

    mytree.Branch('jet4_E',t_var_jet4_E,'jet4_E/F')
    mytree.Branch('jet4_theta',t_var_jet4_theta,'jet4_theta/F')
    mytree.Branch('jet4_BTag',t_var_jet4_BTag,'jet4_BTag/F')
    #mytree.Branch('jet4_jetCharge',t_var_jet4_jetCharge,'jet4_jetCharge/F')

    mytree.Branch('jet5_E',t_var_jet5_E,'jet5_E/F')
    mytree.Branch('jet5_theta',t_var_jet5_theta,'jet5_theta/F')
    mytree.Branch('jet5_BTag',t_var_jet5_BTag,'jet5_BTag/F')
    #mytree.Branch('jet5_jetCharge',t_var_jet5_jetCharge,'jet5_jetCharge/F')

    mytree.Branch('jet6_E',t_var_jet6_E,'jet6_E/F')
    mytree.Branch('jet6_theta',t_var_jet6_theta,'jet6_theta/F')
    mytree.Branch('jet6_BTag',t_var_jet6_BTag,'jet6_BTag/F')
    #mytree.Branch('jet6_jetCharge',t_var_jet6_jetCharge,'jet6_jetCharge/F')

    mytree.Branch('y12',t_var_y12,'y12/F')
    mytree.Branch('y23',t_var_y23,'y23/F')
    mytree.Branch('y34',t_var_y34,'y24/F')
    mytree.Branch('y45',t_var_y45,'y45/F')
    mytree.Branch('y56',t_var_y56,'y56/F')

    mytree.Branch('gen_y12',t_var_gen_y12,'gen_y12/F')
    mytree.Branch('gen_y23',t_var_gen_y23,'gen_y23/F')
    mytree.Branch('gen_y34',t_var_gen_y34,'gen_y34/F')
    mytree.Branch('gen_y45',t_var_gen_y45,'gen_y45/F')
    mytree.Branch('gen_y56',t_var_gen_y56,'gen_y56/F')

    mytree.Branch('BTag1',t_var_BTag1,'BTag1/F')
    mytree.Branch('BTag2',t_var_BTag2,'BTag2/F')
    mytree.Branch('BTag3',t_var_BTag3,'BTag3/F')
    mytree.Branch('BTag4',t_var_BTag4,'BTag4/F')

    mytree.Branch('BTag_sum_max2',t_var_BTag_sum_max2,'BTag_sum_max2/F')
    mytree.Branch('BTag_sum_max3',t_var_BTag_sum_max3,'BTag_sum_max3/F')
    mytree.Branch('BTag_sum_max4',t_var_BTag_sum_max4,'BTag_sum_max4/F')
    mytree.Branch('BTag_sum_all',t_var_BTag_sum_all,'BTag_sum_all/F')


    tree = file.Get("showerData")

    #hist = file2.Get("h_runstatistics")

    #if tree.GetEntries()!=hist.GetBinContent(1):
    #    print "tree_entries not hist content/diff", tree.GetEntries(),hist.GetBinContent(1), tree.GetEntries()-hist.GetBinContent(1)

    sqrtS_low = 750.0
    #sqrtS_high = 2500.0
    sqrtS_high = 2500.0
    #number of leptons need to be smaller than this number
    m_cut_nLeptons = 1
 
    weight = xsec*lumi/tree.GetEntries()
    print "tree-entries ",tree.GetEntries(), " weight ",weight, "xsec",xsec,"lumi",lumi,"total event",weight*tree.GetEntries(), 'to', xsec*lumi, 'size of vector',len(hist_vec_reco_1D_parton)

    num_entry=-1
    num_total_exception=0
    num_count=0

    for ientry in tree:
        num_entry+=1
        #if num_entry<1040:
        #    continue;
        
        #print "in entry ",num_entry

        if num_entry%(int(tree.GetEntries()/5.)) == 0:
            print "in entry ",num_entry
 

        t_var_weight[0] =weight
        t_var_sqrtS[0] = -10.
        t_var_sqrtS_gen[0] = -10.
        t_var_sqrtS_jets [0] = -10.
        t_var_sqrtS_genjets[0] = -10.
        t_var_sqrtS_parton[0]= -10.

        t_var_MET[0] = -10.
        t_var_MET_gen[0] = -10.
        t_var_MHT[0] = -10.
        t_var_MHT_gen[0] = -10.

        
        t_var_comb_jet1_mass[0] = -10.
        t_var_comb_jet1_theta[0] = -10.
        t_var_comb_jet1_phi[0] = -10.
        t_var_comb_jet1_E[0] = -10.
        t_var_comb_jet1_BTagMax[0] = -10.
        #angles between original jets combined in this one
        t_var_comb_jet1_dalpha[0] = -10.
        t_var_comb_jet1_dphi[0] = -10.
        t_var_comb_jet1_dtheta[0] = -10.
        #ratio of more energetic input jet for combination to total sum
        t_var_comb_jet1_E1_over_Etot[0] = -10.
        #t_var_comb_jet1_jetChargeMax[0] = -10.
        #t_var_comb_jet1_jetChargeMin[0] = -10.
        #t_var_comb_jet1_E1_jetCharge[0] = -10.
        
        t_var_comb_jet2_mass[0] = -10.
        t_var_comb_jet2_theta[0] = -10.
        t_var_comb_jet2_phi[0] = -10.
        t_var_comb_jet2_E[0] = -10.
        t_var_comb_jet2_BTagMax[0] = -10.
        t_var_comb_jet2_dalpha[0] = -10.
        t_var_comb_jet2_dphi[0] = -10.
        t_var_comb_jet2_dtheta[0] = -10.
        t_var_comb_jet2_E1_over_Etot[0] = -10.
        #t_var_comb_jet2_jetChargeMax[0] = -10.
        #t_var_comb_jet2_jetChargeMin[0] = -10.
        #t_var_comb_jet2_E1_jetCharge[0] = -10.
        
        t_var_comb_jet3_mass[0] = -10.
        t_var_comb_jet3_theta[0] = -10.
        t_var_comb_jet3_phi[0] = -10.
        t_var_comb_jet3_E[0] = -10.
        t_var_comb_jet3_BTagMax[0] = -10.
        t_var_comb_jet3_dalpha[0] = -10.
        t_var_comb_jet3_dphi[0] = -10.
        t_var_comb_jet3_dtheta[0] = -10.
        t_var_comb_jet3_E1_over_Etot[0] = -10.
        #t_var_comb_jet3_jetChargeMax[0] = -10.
        #t_var_comb_jet3_jetChargeMin[0] = -10.
        #t_var_comb_jet3_E1_jetCharge[0] = -10.
        
        t_var_comb_genjet1_mass[0] = -10.
        t_var_comb_genjet1_theta[0] = -10.
        t_var_comb_genjet1_phi[0] = -10.
        t_var_comb_genjet1_E[0] = -10.
        #angles between original jets combined in this one
        t_var_comb_genjet1_dalpha[0] = -10.
        t_var_comb_genjet1_dphi[0] = -10.
        t_var_comb_genjet1_dtheta[0] = -10.
        #ratio of more energetic input jet for combination to total sum
        t_var_comb_genjet1_E1_over_Etot[0] = -10.
        
        t_var_comb_genjet2_mass[0] = -10.
        t_var_comb_genjet2_theta[0] = -10.
        t_var_comb_genjet2_phi[0] = -10.
        t_var_comb_genjet2_E[0] = -10.
        #angles between original jets combined in this one
        t_var_comb_genjet2_dalpha[0] = -10.
        t_var_comb_genjet2_dphi[0] = -10.
        t_var_comb_genjet2_dtheta[0] = -10.
        #ratio of more energetic input jet for combination to total sum
        t_var_comb_genjet2_E1_over_Etot[0] = -10.
        
        t_var_comb_genjet3_mass[0] = -10.
        t_var_comb_genjet3_theta[0] = -10.
        t_var_comb_genjet3_phi[0] = -10.
        t_var_comb_genjet3_E[0] = -10.
        #angles between original jets combined in this one
        t_var_comb_genjet3_dalpha[0] = -10.
        t_var_comb_genjet3_dphi[0] = -10.
        t_var_comb_genjet3_dtheta[0] = -10.
        #ratio of more energetic input jet for combination to total sum
        t_var_comb_genjet3_E1_over_Etot[0] = -10.
        
        #order jets here by energy
        t_var_jet1_E[0] = -10.
        t_var_jet1_theta[0] = -10.
        t_var_jet1_BTag[0] = -10.
        #t_var_jet1_jetCharge[0] = -10.
        
        t_var_jet2_E[0] = -10.
        t_var_jet2_theta[0] = -10.
        t_var_jet2_BTag[0] = -10.
        #t_var_jet2_jetCharge[0] = -10.
        
        t_var_jet3_E[0] = -10.
        t_var_jet3_theta[0] = -10.
        t_var_jet3_BTag[0] = -10.
        #t_var_jet3_jetCharge[0] = -10.
        
        t_var_jet4_E[0] = -10.
        t_var_jet4_theta[0] = -10.
        t_var_jet4_BTag[0] = -10.
        #t_var_jet4_jetCharge[0] = -10.
        
        t_var_jet5_E[0] = -10.
        t_var_jet5_theta[0] = -10.
        t_var_jet5_BTag[0] = -10.
        #t_var_jet5_jetCharge[0] = -10.
        
        t_var_jet6_E[0] = -10.
        t_var_jet6_theta[0] = -10.
        t_var_jet6_BTag[0] = -10.
        #t_var_jet6_jetCharge[0] = -10.
        
        t_var_y12[0] = -10.
        t_var_y23[0] = -10.
        t_var_y34[0] = -10.
        t_var_y45[0] = -10.
        t_var_y56[0] = -10.
        
        t_var_gen_y12[0] = -10.
        t_var_gen_y23[0] = -10.
        t_var_gen_y34[0] = -10.
        t_var_gen_y45[0] = -10.
        t_var_gen_y56[0] = -10.
        
        #order by BTag --> largest BTag is called BTag1
        t_var_BTag1[0] = -10.
        t_var_BTag2[0] = -10.
        t_var_BTag3[0] = -10.
        t_var_BTag4[0] = -10.
        
        t_var_BTag_sum_max2[0] = -10.
        t_var_BTag_sum_max3[0] = -10.
        t_var_BTag_sum_max4[0] = -10.
        t_var_BTag_sum_all[0] = -10.

        tempTotEventP4=TLorentzVector(0,0,0,0)
        tempTotEventP4HHZ=TLorentzVector(0,0,0,0)
        tempZP4=TLorentzVector(0,0,0,0)
        tempH1P4=TLorentzVector(0,0,0,0)
        tempH2P4=TLorentzVector(0,0,0,0)
        tempHHP4=TLorentzVector(0,0,0,0)
        
        quark_Z_decays=0
        
        tempH1_b=TLorentzVector(0,0,0,0);
        tempH1_bbar=TLorentzVector(0,0,0,0);
        tempH2_b=TLorentzVector(0,0,0,0);
        tempH2_bbar=TLorentzVector(0,0,0,0);

        tempZ_q_pos=TLorentzVector(0,0,0,0)
        tempZ_q_neg=TLorentzVector(0,0,0,0)

        #quarks ordered by energy, no cut on flavour 
        tempH1_q1=TLorentzVector(0,0,0,0)
        tempH1_q2=TLorentzVector(0,0,0,0)

        tempH2_q1=TLorentzVector(0,0,0,0)
        tempH2_q2=TLorentzVector(0,0,0,0)

        tempH1_V1_q1=TLorentzVector(0,0,0,0)
        tempH1_V1_q2=TLorentzVector(0,0,0,0)

        tempH1_V2_q1=TLorentzVector(0,0,0,0)
        tempH1_V2_q2=TLorentzVector(0,0,0,0)

        tempH2_V1_q1=TLorentzVector(0,0,0,0)
        tempH2_V1_q2=TLorentzVector(0,0,0,0)

        tempH2_V2_q1=TLorentzVector(0,0,0,0)
        tempH2_V2_q2=TLorentzVector(0,0,0,0)

        temp_true_inv=TLorentzVector(0,0,0,0)

        true_inv_E=ientry.true_inv_E
        true_inv_Px=ientry.true_inv_Px
        true_inv_Py=ientry.true_inv_Py
        true_inv_Pz=ientry.true_inv_Pz

        temp_true_inv.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E)

        H1H2_decays_bbar = True
        H1H2_decays_qqbar = True
        Z_decays_qqbar = True
        if len(hist_vec_reco_1D_parton)>0:
            trueME_E=ientry.trueME_E
            trueME_Px=ientry.trueME_Px
            trueME_Py=ientry.trueME_Py
            trueME_Pz=ientry.trueME_Pz
            trueME_PDGID=ientry.trueME_PDGID
            trueME_ParentID=ientry.trueME_ParentID
            trueME_PDGID_D1=ientry.trueME_PDGID_D1
            
            index_firstH=-1
            index_secondH=-1
            index_firstZ=-1
        
        #indices 0-1 outgoing electrons after ISR and beam strahlung (i.e. the "real collision")
        #index 2,3 ISR photons of incoming electrons
        #index 4 is Z, index 5 is H1, index 6 is H2
        #index 7 and 8 are daughters from H1
        #index 9 and 10 are decay products of H2
        #index 11 and 12 are decay products of Z 
            for i in range(len(trueME_E)):
                #print 'in ME ',i,trueME_PDGID[i],trueME_E[i],num_entry
                if trueME_PDGID[4]!=25 or trueME_PDGID[5]!=25 or trueME_PDGID[6]!=23:
                    num_total_exception+=1
                if i<2:
                    temp=TLorentzVector(0,0,0,0)
                    temp.SetPxPyPzE(trueME_Px[i],trueME_Py[i],trueME_Pz[i],trueME_E[i])
                    tempTotEventP4+=temp
                if trueME_PDGID[i]==25 and index_firstH==-1:
                    index_firstH=i
                if trueME_PDGID[i]==25 and i!=index_firstH and index_secondH==-1 :
                    index_secondH=i
                if trueME_PDGID[i]==23 and index_firstZ==-1:
                    index_firstZ=i
                #if index_firstZ!=-1 and index_firstH!=-1 and index_secondH!=-1 :
                #    break


            index_firstH_d1=index_firstH+2
            index_firstH_d2=index_firstH+3

            index_secondH_d1=index_secondH+3
            index_secondH_d2=index_secondH+4

            indexZ_d1=index_firstZ+7
            indexZ_d2=index_firstZ+8

            if ishhzfile:
                #here the order is Z, then H and H
                index_firstH_d1=index_firstH+3
                index_firstH_d2=index_firstH+4
                index_secondH_d1=index_secondH+4
                index_secondH_d2=index_secondH+5                
                indexZ_d1=index_firstZ+5
                indexZ_d2=index_firstZ+6

            if (abs(trueME_PDGID[index_firstH_d1])!=5 or abs(trueME_PDGID[index_firstH_d2])!=5) or (abs(trueME_PDGID[index_secondH_d1])!=5 or abs(trueME_PDGID[index_secondH_d1])!=5):
                H1H2_decays_bbar=False

            if (abs(trueME_PDGID[index_firstH_d1])>6 or abs(trueME_PDGID[index_firstH_d2])>6) or (abs(trueME_PDGID[index_secondH_d1])>6 or abs(trueME_PDGID[index_secondH_d1])>6):
                H1H2_decays_qqbar=False

        #H1 decays into quarks or gluons, aka jets
            if abs(trueME_PDGID[index_firstH_d1])<7 or trueME_PDGID[index_firstH_d1]==21 :
                if trueME_ParentID[index_firstH_d1]!=1 or trueME_ParentID[index_firstH_d2]!=1:
                    print 'should have been daughters from H1',trueME_ParentID[index_firstH_d1],trueME_ParentID[index_firstH_d2]
                if trueME_E[index_firstH_d1]> trueME_E[index_firstH_d2]:
                    tempH1_q1.SetPxPyPzE(trueME_Px[index_firstH_d1],trueME_Py[index_firstH_d1],trueME_Pz[index_firstH_d1],trueME_E[index_firstH_d1])
                    tempH1_q2.SetPxPyPzE(trueME_Px[index_firstH_d2],trueME_Py[index_firstH_d2],trueME_Pz[index_firstH_d2],trueME_E[index_firstH_d2])
                else :
                    tempH1_q2.SetPxPyPzE(trueME_Px[index_firstH_d1],trueME_Py[index_firstH_d1],trueME_Pz[index_firstH_d1],trueME_E[index_firstH_d1])
                    tempH1_q1.SetPxPyPzE(trueME_Px[index_firstH_d2],trueME_Py[index_firstH_d2],trueME_Pz[index_firstH_d2],trueME_E[index_firstH_d2])
                if trueME_PDGID[index_firstH_d1]==5 :
                    tempH1_b.SetPxPyPzE(trueME_Px[index_firstH_d1],trueME_Py[index_firstH_d1],trueME_Pz[index_firstH_d1],trueME_E[index_firstH_d1])
                    tempH1_bbar.SetPxPyPzE(trueME_Px[index_firstH_d2],trueME_Py[index_firstH_d2],trueME_Pz[index_firstH_d2],trueME_E[index_firstH_d2])
                    if abs(trueME_PDGID[index_firstH_d1]) != abs(trueME_PDGID[index_firstH_d2]):
                        print 'pdg id should for H1 be the same in abs value, what is wrong ',trueME_PDGID[index_firstH_d1],-trueME_PDGID[index_firstH_d2]
                elif trueME_PDGID[index_firstH_d1]==-5 :
                    tempH1_bbar.SetPxPyPzE(trueME_Px[index_firstH_d1],trueME_Py[index_firstH_d1],trueME_Pz[index_firstH_d1],trueME_E[index_firstH_d1])
                    tempH1_b.SetPxPyPzE(trueME_Px[index_firstH_d2],trueME_Py[index_firstH_d2],trueME_Pz[index_firstH_d2],trueME_E[index_firstH_d2])             
                if trueME_PDGID[index_firstH_d1] != -trueME_PDGID[index_firstH_d2]:
                    if H1H2_decays_bbar:
                        print 'pdg id should for H1 be the same, what is wrong bbbar',trueME_PDGID[index_firstH_d1],-trueME_PDGID[index_firstH_d2]
        #H2 decays into quarks or gluons, aka jets (after index, +1 is Z, +2 H1 d1 +3 --> H1 d2
            if abs(trueME_PDGID[index_secondH_d1])<7 or trueME_PDGID[index_secondH_d2]==21 :
                if trueME_ParentID[index_secondH_d1]!=2 or trueME_ParentID[index_secondH_d2]!=2:
                    print 'should have been daughters from H2',trueME_ParentID[index_firstH_d1],trueME_ParentID[index_firstH_d2],num_entry
                if trueME_E[index_secondH_d1]> trueME_E[index_secondH_d2]:
                    tempH2_q1.SetPxPyPzE(trueME_Px[index_secondH_d1],trueME_Py[index_secondH_d1],trueME_Pz[index_secondH_d1],trueME_E[index_secondH_d1])
                    tempH2_q2.SetPxPyPzE(trueME_Px[index_secondH_d2],trueME_Py[index_secondH_d2],trueME_Pz[index_secondH_d2],trueME_E[index_secondH_d2])
                else :
                    tempH2_q2.SetPxPyPzE(trueME_Px[index_secondH_d1],trueME_Py[index_secondH_d1],trueME_Pz[index_secondH_d1],trueME_E[index_secondH_d1])
                    tempH2_q1.SetPxPyPzE(trueME_Px[index_secondH_d2],trueME_Py[index_secondH_d2],trueME_Pz[index_secondH_d2],trueME_E[index_secondH_d2])
                if trueME_PDGID[index_secondH_d1]==5 :
                    tempH2_b.SetPxPyPzE(trueME_Px[index_secondH_d1],trueME_Py[index_secondH_d1],trueME_Pz[index_secondH_d1],trueME_E[index_secondH_d1])
                    tempH2_bbar.SetPxPyPzE(trueME_Px[index_secondH_d2],trueME_Py[index_secondH_d2],trueME_Pz[index_secondH_d2],trueME_E[index_secondH_d2])
                    if abs(trueME_PDGID[index_secondH_d1]) != abs(trueME_PDGID[index_secondH_d2]):
                        print 'pdg id should for H2 be the same, what is wrong ',trueME_PDGID[index_secondH_d1],-trueME_PDGID[index_secondH_d2]
                    elif trueME_PDGID[index_secondH_d1]==-5 :
                        tempH2_bbar.SetPxPyPzE(trueME_Px[index_secondH_d1],trueME_Py[index_secondH_d1],trueME_Pz[index_secondH_d1],trueME_E[index_secondH_d1])
                        tempH2_b.SetPxPyPzE(trueME_Px[index_secondH_d2],trueME_Py[index_secondH_d2],trueME_Pz[index_secondH_d2],trueME_E[index_secondH_d2])             
                if trueME_PDGID[index_secondH_d1] != -trueME_PDGID[index_secondH_d2]:
                    if H1H2_decays_bbar:
                        print 'pdg id should for H2 be the same, what is wrong bbbar ',trueME_PDGID[index_secondH_d1],-trueME_PDGID[index_secondH_d2]

            if (trueME_PDGID[indexZ_d1]>6 or  trueME_PDGID[indexZ_d2]>6):
            #print 'index of Z daughter is',trueME_PDGID[indexZ_d1]
                Z_decays_qqbar=False
        
        #positive charge, aka up type quarks or down bar type quarks
            if(trueME_PDGID[indexZ_d1]==-1 or trueME_PDGID[indexZ_d1]==2 or trueME_PDGID[indexZ_d1]==-3 or trueME_PDGID[indexZ_d1]==4 or trueME_PDGID[indexZ_d1]==-5 or trueME_PDGID[indexZ_d1]==6):
                if trueME_PDGID[indexZ_d1] != -trueME_PDGID[indexZ_d2]:
                    print 'pdg id should be the same for Z, what is wrong ',trueME_PDGID[indexZ_d1],-trueME_PDGID[indexZ_d2]
                if trueME_ParentID[indexZ_d1]!=3 or trueME_ParentID[indexZ_d2]!=3:
                    print 'should have been daughters from Z',trueME_ParentID[indexZ_d1],trueME_ParentID[indexZ_d2]
                tempZ_q_pos.SetPxPyPzE(trueME_Px[indexZ_d1],trueME_Py[indexZ_d1],trueME_Pz[indexZ_d1],trueME_E[indexZ_d1])
                tempZ_q_neg.SetPxPyPzE(trueME_Px[indexZ_d2],trueME_Py[indexZ_d2],trueME_Pz[indexZ_d2],trueME_E[indexZ_d2])
            else:
            #quark index 6 is negatively charged
                if trueME_PDGID[indexZ_d1] != -trueME_PDGID[indexZ_d2]:
                    print 'pdg id should be the same for Z, what is wrong ',trueME_PDGID[indexZ_d1],-trueME_PDGID[indexZ_d2]
                if trueME_ParentID[indexZ_d1]!=3 or trueME_ParentID[indexZ_d2]!=3:
                    print 'should have been daughters from Z',trueME_ParentID[indexZ_d1],trueME_ParentID[indexZ_d2]
                tempZ_q_neg.SetPxPyPzE(trueME_Px[indexZ_d1],trueME_Py[indexZ_d1],trueME_Pz[indexZ_d1],trueME_E[indexZ_d1])
                tempZ_q_pos.SetPxPyPzE(trueME_Px[indexZ_d2],trueME_Py[indexZ_d2],trueME_Pz[indexZ_d2],trueME_E[indexZ_d2])
    


 
            
            tempZ_first=TLorentzVector(0,0,0,0)
            tempZ_first.SetPxPyPzE(trueME_Px[index_firstZ],trueME_Py[index_firstZ],trueME_Pz[index_firstZ],trueME_E[index_firstZ])

            if(trueME_E[index_firstH]>trueME_E[index_secondH]):
                tempH1P4.SetPxPyPzE(trueME_Px[index_firstH],trueME_Py[index_firstH],trueME_Pz[index_firstH],trueME_E[index_firstH])
                tempH2P4.SetPxPyPzE(trueME_Px[index_secondH],trueME_Py[index_secondH],trueME_Pz[index_secondH],trueME_E[index_secondH])
            else: 
                tempH2P4.SetPxPyPzE(trueME_Px[index_firstH],trueME_Py[index_firstH],trueME_Pz[index_firstH],trueME_E[index_firstH])
                tempH1P4.SetPxPyPzE(trueME_Px[index_secondH],trueME_Py[index_secondH],trueME_Pz[index_secondH],trueME_E[index_secondH])

            tempHHP4=tempH1P4+tempH2P4
            tempZP4=tempZ_q_pos+tempZ_q_neg
            tempTotEventP4HHZ=tempHHP4+tempZ_first
        #print "in entry ",num_entry,tempTotEventP4.M()
            hist_vec_reco_1D_parton[0].Fill(tempTotEventP4.M(),weight)
            hist_vec_reco_1D_parton[42].Fill(tempTotEventP4HHZ.M(),weight)
            hist_vec_reco_1D_parton[1].Fill(tempH1P4.E()+tempH2P4.E(),weight)

            if(tempH1P4.E()>tempH2P4.E()):
                hist_vec_reco_1D_parton[2].Fill(tempH1P4.E(),weight)
                hist_vec_reco_1D_parton[3].Fill(tempH2P4.E(),weight)
                hist_vec_reco_1D_parton[118].Fill(tempH1P4.P(),weight)
                hist_vec_reco_1D_parton[119].Fill(tempH2P4.P(),weight)
                hist_vec_reco_1D_parton[62].Fill(degrees(tempH1P4.Theta()),weight)
                hist_vec_reco_1D_parton[63].Fill(degrees(tempH2P4.Theta()),weight)
                hist_vec_reco_1D_parton[121].Fill(1-cos(tempH1P4.Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[122].Fill(1-cos(tempH2P4.Angle(tempZP4.Vect())),weight)
            else:
                print 'we should NOT have this case anymore'
                hist_vec_reco_1D_parton[2].Fill(tempH2P4.E(),weight)
                hist_vec_reco_1D_parton[3].Fill(tempH1P4.E(),weight)
                hist_vec_reco_1D_parton[118].Fill(tempH2P4.P(),weight)
                hist_vec_reco_1D_parton[119].Fill(tempH1P4.P(),weight)
                hist_vec_reco_1D_parton[62].Fill(degrees(tempH2P4.Theta()),weight)
                hist_vec_reco_1D_parton[63].Fill(degrees(tempH1P4.Theta()),weight)
                hist_vec_reco_1D_parton[121].Fill(1-cos(tempH2P4.Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[122].Fill(1-cos(tempH1P4.Angle(tempZP4.Vect())),weight)
            hist_vec_reco_1D_parton[4].Fill(tempZ_first.E(),weight)
            hist_vec_reco_1D_parton[120].Fill(tempZ_first.P(),weight)

            hist_vec_reco_1D_parton[5].Fill((tempH1P4+tempH2P4).M(),weight)
            hist_vec_reco_1D_parton[6].Fill(degrees((tempH1P4+tempH2P4).Angle(tempZ_first.Vect())),weight)
            hist_vec_reco_1D_parton[7].Fill(degrees(DeltaPhi((tempH1P4+tempH2P4).Phi(),tempZ_first.Phi())),weight)
            hist_vec_reco_1D_parton[8].Fill(degrees(abs((tempH1P4+tempH2P4).Theta()-tempZ_first.Theta())),weight)
            hist_vec_reco_1D_parton[9].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D_parton[10].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
            hist_vec_reco_1D_parton[11].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
            hist_vec_reco_1D_parton[39].Fill(1.-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D_parton[64].Fill(degrees(tempZP4.Theta()),weight)

        #for histogram filling: H1 energy> H2 energy
            if H1H2_decays_bbar:
                if(tempH1P4.E()>tempH2P4.E()):
                    hist_vec_reco_1D_parton[12].Fill(degrees(tempH1_bbar.Angle(tempH1_b.Vect())),weight)
                    hist_vec_reco_1D_parton[13].Fill(degrees(DeltaPhi((tempH1_bbar).Phi(),tempH1_b.Phi())),weight)
                    hist_vec_reco_1D_parton[14].Fill(degrees(abs((tempH1_bbar).Theta()-tempH1_b.Theta())),weight)
                    hist_vec_reco_1D_parton[60].Fill(1.-cos(tempH1_bbar.Angle(tempH1_b.Vect())),weight)

                    hist_vec_reco_1D_parton[15].Fill(degrees(tempH2_bbar.Angle(tempH2_b.Vect())),weight)
                    hist_vec_reco_1D_parton[16].Fill(degrees(DeltaPhi((tempH2_bbar).Phi(),tempH2_b.Phi())),weight)
                    hist_vec_reco_1D_parton[17].Fill(degrees(abs((tempH2_bbar).Theta()-tempH2_b.Theta())),weight)
                    hist_vec_reco_1D_parton[61].Fill(1.-cos(tempH2_bbar.Angle(tempH2_b.Vect())),weight)
                else:
                    hist_vec_reco_1D_parton[12].Fill(degrees(tempH2_bbar.Angle(tempH2_b.Vect())),weight)
                    hist_vec_reco_1D_parton[13].Fill(degrees(DeltaPhi((tempH2_bbar).Phi(),tempH2_b.Phi())),weight)
                    hist_vec_reco_1D_parton[14].Fill(degrees(abs(tempH2_bbar.Theta()-tempH2_b.Theta())),weight)
                    hist_vec_reco_1D_parton[60].Fill(1.-cos(tempH2_bbar.Angle(tempH2_b.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[15].Fill(degrees(tempH1_bbar.Angle(tempH1_b.Vect())),weight)
                    hist_vec_reco_1D_parton[16].Fill(degrees(DeltaPhi((tempH1_bbar).Phi(),tempH1_b.Phi())),weight)
                    hist_vec_reco_1D_parton[17].Fill(degrees(abs((tempH1_bbar).Theta()-tempH1_b.Theta())),weight)
                    hist_vec_reco_1D_parton[61].Fill(1.-cos(tempH1_bbar.Angle(tempH1_b.Vect())),weight)


            if H1H2_decays_qqbar:
                if(tempH1P4.E()>tempH2P4.E()):
                    hist_vec_reco_1D_parton[18].Fill(degrees(tempH1_q2.Angle(tempH1_q1.Vect())),weight)
                    hist_vec_reco_1D_parton[19].Fill(degrees(DeltaPhi((tempH1_q2).Phi(),tempH1_q1.Phi())),weight)
                    hist_vec_reco_1D_parton[20].Fill(degrees(abs((tempH1_q2).Theta()-tempH1_q1.Theta())),weight)
                    hist_vec_reco_1D_parton[40].Fill(1.-cos(tempH1_q2.Angle(tempH1_q1.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[21].Fill(degrees(tempH2_q2.Angle(tempH2_q1.Vect())),weight)
                    hist_vec_reco_1D_parton[22].Fill(degrees(DeltaPhi((tempH2_q2).Phi(),tempH2_q1.Phi())),weight)
                    hist_vec_reco_1D_parton[23].Fill(degrees(abs((tempH2_q2).Theta()-tempH2_q1.Theta())),weight)
                    hist_vec_reco_1D_parton[41].Fill(1.-cos(tempH2_q2.Angle(tempH2_q1.Vect())),weight)
                else:
                    hist_vec_reco_1D_parton[18].Fill(degrees(tempH2_q2.Angle(tempH2_q1.Vect())),weight)
                    hist_vec_reco_1D_parton[19].Fill(degrees(DeltaPhi((tempH2_q2).Phi(),tempH2_q1.Phi())),weight)
                    hist_vec_reco_1D_parton[20].Fill(degrees(abs(tempH2_q2.Theta()-tempH2_q1.Theta())),weight)
                    hist_vec_reco_1D_parton[40].Fill(1.-cos(tempH2_q2.Angle(tempH2_q1.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[21].Fill(degrees(tempH1_q2.Angle(tempH1_q1.Vect())),weight)
                    hist_vec_reco_1D_parton[22].Fill(degrees(DeltaPhi((tempH1_q2).Phi(),tempH1_q1.Phi())),weight)
                    hist_vec_reco_1D_parton[23].Fill(degrees(abs((tempH1_q2).Theta()-tempH1_q1.Theta())),weight)
                    hist_vec_reco_1D_parton[41].Fill(1.-cos(tempH1_q2.Angle(tempH1_q1.Vect())),weight)

                if Z_decays_qqbar:
                    hist_vec_reco_1D_parton[24].Fill(degrees(tempZ_q_neg.Angle(tempZ_q_pos.Vect())),weight)
                    hist_vec_reco_1D_parton[25].Fill(degrees(DeltaPhi((tempZ_q_neg).Phi(),tempZ_q_pos.Phi())),weight)
                    hist_vec_reco_1D_parton[26].Fill(degrees(abs((tempZ_q_neg).Theta()-tempZ_q_pos.Theta())),weight)
                    hist_vec_reco_1D_parton[59].Fill(1.-cos(tempZ_q_neg.Angle(tempZ_q_pos.Vect())),weight)


            tempHClosest=TLorentzVector(0,0,0,0)
            tempHClosest=tempH1P4
            tempHFar=TLorentzVector(0,0,0,0)
            tempHFar=tempH2P4
            if tempH2P4.Angle(tempZ_first.Vect())<tempH1P4.Angle(tempZ_first.Vect()):
                hist_vec_reco_1D_parton[43].Fill(degrees(tempZ_first.Angle(tempH2P4.Vect())),weight)
                hist_vec_reco_1D_parton[44].Fill(degrees(DeltaPhi(tempZ_first.Phi(),tempH2P4.Phi())),weight)
                hist_vec_reco_1D_parton[45].Fill(degrees(abs(tempZ_first.Theta()-tempH2P4.Theta())),weight)
                hist_vec_reco_1D_parton[46].Fill(degrees(tempZ_first.Angle(tempH1P4.Vect())),weight)
                hist_vec_reco_1D_parton[47].Fill(degrees(DeltaPhi(tempZ_first.Phi(),tempH1P4.Phi())),weight)
                hist_vec_reco_1D_parton[48].Fill(degrees(abs(tempZ_first.Theta()-tempH1P4.Theta())),weight)
                hist_vec_reco_1D_parton[49].Fill((tempZ_first+tempH2P4).E(),weight)
                hist_vec_reco_1D_parton[50].Fill((tempZ_first+tempH2P4).E()-tempH1P4.E(),weight)
                
                tempHClosest=tempH2P4
                tempHFar=tempH1P4
            #checked that H2 closer to Z than H1 to Z
                if tempH2P4.Angle(tempH1P4.Vect())< tempH2P4.Angle(tempZ_first.Vect()):
                #H1 closer to H2 than H2 to Z
                    hist_vec_reco_1D_parton[51].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())),weight)
                    hist_vec_reco_1D_parton[52].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
                    hist_vec_reco_1D_parton[53].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
                    hist_vec_reco_1D_parton[57].Fill(1-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
                else:
                #H2 to Z closest angle
                    hist_vec_reco_1D_parton[51].Fill(degrees(tempH2P4.Angle(tempZ_first.Vect())),weight)
                    hist_vec_reco_1D_parton[52].Fill(degrees(DeltaPhi(tempH2P4.Phi(),tempZ_first.Phi())),weight)
                    hist_vec_reco_1D_parton[53].Fill(degrees(abs(tempH2P4.Theta()-tempZ_first.Theta())),weight)
                    hist_vec_reco_1D_parton[57].Fill(1-cos(tempZ_first.Angle(tempH2P4.Vect())),weight)
                #H1 and H2 further appart than H1 to Z
                    if tempH2P4.Angle(tempH1P4.Vect())> tempH1P4.Angle(tempZ_first.Vect()):
                    #H1 and H2 further appart than H1 to Z
                        hist_vec_reco_1D_parton[54].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())),weight)
                        hist_vec_reco_1D_parton[55].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
                        hist_vec_reco_1D_parton[56].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
                        hist_vec_reco_1D_parton[58].Fill(1-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
                    else:
                    #then H1 to Z > H1,H2 > H2,Z
                        hist_vec_reco_1D_parton[54].Fill(degrees(tempH1P4.Angle(tempZ_first.Vect())))
                        hist_vec_reco_1D_parton[55].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempZ_first.Phi())),weight)
                        hist_vec_reco_1D_parton[56].Fill(degrees(abs(tempH1P4.Theta()-tempZ_first.Theta())),weight)
                        hist_vec_reco_1D_parton[58].Fill(1-cos(tempH1P4.Angle(tempZ_first.Vect())),weight)
            else:
            #clear H1 closer to Z than H2 to Z
                hist_vec_reco_1D_parton[43].Fill(degrees(tempZ_first.Angle(tempH1P4.Vect())),weight)
                hist_vec_reco_1D_parton[44].Fill(degrees(DeltaPhi(tempZ_first.Phi(),tempH1P4.Phi())),weight)
                hist_vec_reco_1D_parton[45].Fill(degrees(abs(tempZ_first.Theta()-tempH1P4.Theta())),weight)
                hist_vec_reco_1D_parton[46].Fill(degrees(tempZ_first.Angle(tempH2P4.Vect())),weight)
                hist_vec_reco_1D_parton[47].Fill(degrees(DeltaPhi(tempZ_first.Phi(),tempH2P4.Phi())),weight)
                hist_vec_reco_1D_parton[48].Fill(degrees(abs(tempZ_first.Theta()-tempH2P4.Theta())),weight)
                hist_vec_reco_1D_parton[49].Fill((tempZ_first+tempH1P4).E(),weight)
                hist_vec_reco_1D_parton[50].Fill((tempZ_first+tempH1P4).E()-tempH2P4.E(),weight)
                if tempH2P4.Angle(tempH1P4.Vect())< tempH1P4.Angle(tempZ_first.Vect()):
                #H1 closer to H2 than H1 to Z
                    hist_vec_reco_1D_parton[51].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())))
                    hist_vec_reco_1D_parton[52].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
                    hist_vec_reco_1D_parton[53].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
                    hist_vec_reco_1D_parton[57].Fill(1-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
                else:
                #H1 to Z closest angle
                    hist_vec_reco_1D_parton[51].Fill(degrees(tempH1P4.Angle(tempZ_first.Vect())))
                    hist_vec_reco_1D_parton[52].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempZ_first.Phi())),weight)
                    hist_vec_reco_1D_parton[53].Fill(degrees(abs(tempH1P4.Theta()-tempZ_first.Theta())),weight)
                    hist_vec_reco_1D_parton[57].Fill(1-cos(tempH1P4.Angle(tempZ_first.Vect())),weight)
                #H1 and H2 further appart than H2 to Z
                    if tempH2P4.Angle(tempH1P4.Vect())> tempH2P4.Angle(tempZ_first.Vect()):
                    #H1 and H2 further appart than H2 to Z, and H2,Z > H1,Z
                        hist_vec_reco_1D_parton[54].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())))
                        hist_vec_reco_1D_parton[55].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
                        hist_vec_reco_1D_parton[56].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
                        hist_vec_reco_1D_parton[58].Fill(1-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
                    else:
                    #then H2 to Z > H1,H2 > H1,Z
                        hist_vec_reco_1D_parton[54].Fill(degrees(tempH2P4.Angle(tempZ_first.Vect())))
                        hist_vec_reco_1D_parton[55].Fill(degrees(DeltaPhi(tempH2P4.Phi(),tempZ_first.Phi())),weight)
                        hist_vec_reco_1D_parton[56].Fill(degrees(abs(tempH2P4.Theta()-tempZ_first.Theta())),weight)
                        hist_vec_reco_1D_parton[58].Fill(1-cos(tempH2P4.Angle(tempZ_first.Vect())),weight)


        #fill now bbar only histos
        if(H1H2_decays_bbar==False):
            continue;
        if(Z_decays_qqbar==False):
            continue;
        num_count+=1
        #print 'H/H/Z E',tempH1P4.E(),tempH2P4.E(),tempZP4.E()
        #print 'H/H/Z P',tempH1P4.P(),tempH2P4.P(),tempZP4.P()
            
        genjet_E=ientry.genjet_E
        genjet_Px=ientry.genjet_Px
        genjet_Py=ientry.genjet_Py
        genjet_Pz=ientry.genjet_Pz
        
        temp_gj_sum=TLorentzVector(0,0,0,0)
        
        for ind in range(len(genjet_E)):
            temp_gj=TLorentzVector(0,0,0,0)
            temp_gj.SetPxPyPzE(genjet_Px[ind],genjet_Py[ind],genjet_Pz[ind],genjet_E[ind])
            temp_gj_sum+=temp_gj
            #print' inex/E/sumE ',ind,temp_gj.E(),temp_gj_sum.E()

            
        recojet_E=ientry.recojet_E
        recojet_Px=ientry.recojet_Px
        recojet_Py=ientry.recojet_Py
        recojet_Pz=ientry.recojet_Pz
        
        #if len(recojet_E)<6:
        #    print 'too small recojet length',len(recojetE)<3
       #     continue;
        recojet_vector=[]
            
        for ind in range(len(recojet_E)):
            temp_=TLorentzVector(0,0,0,0)
            temp_.SetPxPyPzE(recojet_Px[ind],recojet_Py[ind],recojet_Pz[ind],recojet_E[ind])
            recojet_vector.append(temp_)
                
        recojet_vector.sort(key=lambda x: x.M(), reverse=True)

        recojet_rfj_E=ientry.recojet_subjet_rfj_j_E
        recojet_rfj_Px=ientry.recojet_subjet_rfj_j_Px
        recojet_rfj_Py=ientry.recojet_subjet_rfj_j_Py
        recojet_rfj_Pz=ientry.recojet_subjet_rfj_j_Pz
        recojet_rfj_BTag=ientry.recojet_subjet_rfj_j_BTag
        #recojet_rfj_jetcharge=ientry.recojet_subjet_jetChargeE_kappa_0_30

        if len(recojet_rfj_E)<len(recojet_E):
            print 'rfjet and jet size were supposed to be the same',len(recojet_rfj_E),len(recojet_E)
        if len(recojet_rfj_E)<6:
            continue
        #if len(recojet_rfj_jetcharge)<len(recojet_rfj_E):
        #    print 'rfjet and jet size were supposed to be the same jet charge vs rfjet',len(recojet_rfj_E),len(recojet_rfj_jetcharge)
            #continue

        genjet_vector=[]
        tot_genjet_Vector_=TLorentzVector(0,0,0,0)
        if len(genjet_E)<3:
            print 'too small genjet length',len(genjetE)<3            
        for ind in range(len(genjet_E)):
            temp_=TLorentzVector(0,0,0,0)
            temp_.SetPxPyPzE(genjet_Px[ind],genjet_Py[ind],genjet_Pz[ind],genjet_E[ind])
            genjet_vector.append(temp_)
            tot_genjet_Vector_+=temp_
        #for vect in range(len(genjet_vector)):
        #    print 'gen mass vect before ',vect,genjet_vector[vect].M()
        
        genjet_vector.sort(key=lambda x: x.M(), reverse=True)
        #for vect in range(len(genjet_vector)):
        #    print num_entry,'gen mass vect after, index ',vect,genjet_vector[vect].M(),genjet_vector[vect].E()
        if len(genjet_E)>3:
            genjet_vector_combto3,ind_gj1,ind_gj2=map(list,zip(*orderedThreeVector(genjet_vector)))
        #for vect in genjet_vector_combto3:
        #   print num_entry,'gen mass vect after combination to 3 ',vect.M(),vect.E()

        recojet_rfj_BTag_Set=set()
        recojet_rfj_vector=[]
        index_rfj_match=[]
        index_rfj_match_set=set()
        tot_rfjet_Vector_=TLorentzVector(0,0,0,0)
        #put the BTag information from rfjet into a combined set with lorentzvector based on rfjet four vector information
        for ind in range(len(recojet_rfj_E)):
            temp_rfj_=TLorentzVector(0,0,0,0)
            temp_rfj_.SetPxPyPzE(recojet_rfj_Px[ind],recojet_rfj_Py[ind],recojet_rfj_Pz[ind],recojet_rfj_E[ind])
            recojet_rfj_BTag_Set.add((temp_rfj_,recojet_rfj_BTag[ind]))
            #print ' in rfjet filling',ind,recojet_rfj_E[ind],recojet_rfj_BTag[ind]
            recojet_rfj_vector.append(temp_rfj_)
            tot_rfjet_Vector_+=temp_rfj_
            #print 'rfjet E/sumE',ind,temp_rfj_.E(),tot_rfjet_Vector_.E()
            angle_min = 360.
            ind_match=-1
        #try not to match rfjets to recovectors (from fastjet) for later processing
            for ind2 in range(len(recojet_vector)):
                angle_temp=degrees(temp_rfj_.Angle(recojet_vector[ind2].Vect()));
                if(angle_temp<angle_min):
                    angle_min=degrees(temp_rfj_.Angle(recojet_vector[ind2].Vect()));
                    ind_match=ind2
            index_rfj_match.append(ind_match)
            index_rfj_match_set.add(ind_match)

        recojet_rfj_BTag_Set_sorted=sorted(recojet_rfj_BTag_Set, key=lambda x:x[1], reverse=True)

 


        #don't use rfjet BTag info here
        #recojet_rfj_vector.sort(key=lambda x: x.M(), reverse=True)
        #for 6 jets order doesn't make sense BUT if we order here, then BTagging info totally screwed up
        recojet_rfj_vector_combto3,ind_rj_rfj_1,ind_rj_rfj_2=map(list,zip(*orderedThreeVector(recojet_rfj_vector)))
        #for vect in range(len(recojet_rfj_vector_combto3)):
        #    print num_entry,'reco _rfj mass vect after combination to 3 ',vect,recojet_rfj_vector_combto3[vect].M(),recojet_rfj_vector_combto3[vect].E(),(recojet_rfj_vector[ind_rj_rfj_1[vect]]+recojet_rfj_vector[ind_rj_rfj_2[vect]]).M(),(recojet_rfj_vector[ind_rj_rfj_1[vect]]+recojet_rfj_vector[ind_rj_rfj_2[vect]]).E()
        #for rfj in range(len(recojet_rfj_vector_combto3)):
        #    print 'rfj M,E, rjf1 E, rfj2 E, ind1/ind2',recojet_rfj_vector_combto3[rfj].M(),recojet_rfj_vector_combto3[rfj].E(),ind_rj_rfj_1[rfj],ind_rj_rfj_2[rfj],len(recojet_rfj_jetcharge)
        #    print 'rfj M,E, rjf1 E, rfj2 E, rfj1 B, rfj2 B',recojet_rfj_vector_combto3[rfj].M(),recojet_rfj_vector_combto3[rfj].E(),recojet_rfj_vector[ind_rj_rfj_1[rfj]].E(),recojet_rfj_BTag[ind_rj_rfj_1[rfj]],recojet_rfj_jetcharge[ind_rj_rfj_1[rfj]],recojet_rfj_vector[ind_rj_rfj_2[rfj]].E(),recojet_rfj_BTag[ind_rj_rfj_2[rfj]],recojet_rfj_jetcharge[ind_rj_rfj_2[rfj]],max(recojet_rfj_jetcharge[ind_rj_rfj_1[rfj]],recojet_rfj_jetcharge[ind_rj_rfj_2[rfj]],key=abs),min(recojet_rfj_jetcharge[ind_rj_rfj_1[rfj]],recojet_rfj_jetcharge[ind_rj_rfj_2[rfj]],key=abs)


        t_var_comb_jet1_mass[0] = recojet_rfj_vector_combto3[0].M()
        t_var_comb_jet1_theta[0] = recojet_rfj_vector_combto3[0].Theta()
        t_var_comb_jet1_phi[0] = recojet_rfj_vector_combto3[0].Phi()
        t_var_comb_jet1_E[0] = recojet_rfj_vector_combto3[0].E()
        t_var_comb_jet1_BTagMax[0] = max(recojet_rfj_BTag[ind_rj_rfj_1[0]],recojet_rfj_BTag[ind_rj_rfj_2[0]])
        #angles between original jets combined in this one
        t_var_comb_jet1_dalpha[0] = recojet_rfj_vector[ind_rj_rfj_1[0]].Angle(recojet_rfj_vector[ind_rj_rfj_2[0]].Vect())
        t_var_comb_jet1_dphi[0] = recojet_rfj_vector[ind_rj_rfj_1[0]].DeltaPhi(recojet_rfj_vector[ind_rj_rfj_2[0]])
        t_var_comb_jet1_dtheta[0] = abs(recojet_rfj_vector[ind_rj_rfj_1[0]].Theta()-recojet_rfj_vector[ind_rj_rfj_2[0]].Theta())
        #ratio of more energetic input jet for combination to total sum
        t_var_comb_jet1_E1_over_Etot[0] = max(recojet_rfj_vector[ind_rj_rfj_1[0]].E(),recojet_rfj_vector[ind_rj_rfj_2[0]].E())/recojet_rfj_vector_combto3[0].E()
        #t_var_comb_jet1_jetChargeMax[0] = max(recojet_rfj_jetcharge[ind_rj_rfj_1[0]],recojet_rfj_jetcharge[ind_rj_rfj_2[0]],key=abs)
        #t_var_comb_jet1_jetChargeMin[0] = min(recojet_rfj_jetcharge[ind_rj_rfj_1[0]],recojet_rfj_jetcharge[ind_rj_rfj_2[0]],key=abs)
        #if recojet_rfj_vector[ind_rj_rfj_1[0]].E()>recojet_rfj_vector[ind_rj_rfj_2[0]].E():
        #    t_var_comb_jet1_E1_jetCharge[0] = recojet_rfj_jetcharge[ind_rj_rfj_1[0]]
        #else:
        #    t_var_comb_jet1_E1_jetCharge[0] = recojet_rfj_jetcharge[ind_rj_rfj_2[0]]

        t_var_comb_jet2_mass[0] = recojet_rfj_vector_combto3[1].M()
        t_var_comb_jet2_theta[0] = recojet_rfj_vector_combto3[1].Theta()
        t_var_comb_jet2_phi[0] = recojet_rfj_vector_combto3[1].Phi()
        t_var_comb_jet2_E[0] = recojet_rfj_vector_combto3[1].E()
        t_var_comb_jet2_BTagMax[0] = max(recojet_rfj_BTag[ind_rj_rfj_1[1]],recojet_rfj_BTag[ind_rj_rfj_2[1]])
        #angles between original jets combined in this one
        t_var_comb_jet2_dalpha[0] = recojet_rfj_vector[ind_rj_rfj_1[1]].Angle(recojet_rfj_vector[ind_rj_rfj_2[1]].Vect())
        t_var_comb_jet2_dphi[0] = recojet_rfj_vector[ind_rj_rfj_1[1]].DeltaPhi(recojet_rfj_vector[ind_rj_rfj_2[1]])
        t_var_comb_jet2_dtheta[0] = abs(recojet_rfj_vector[ind_rj_rfj_1[1]].Theta()-recojet_rfj_vector[ind_rj_rfj_2[1]].Theta())
        #ratio of more energetic input jet for combination to total sum
        t_var_comb_jet2_E1_over_Etot[0] = max(recojet_rfj_vector[ind_rj_rfj_1[1]].E(),recojet_rfj_vector[ind_rj_rfj_2[1]].E())/recojet_rfj_vector_combto3[1].E()
        #t_var_comb_jet2_jetChargeMax[0] = max(recojet_rfj_jetcharge[ind_rj_rfj_1[1]],recojet_rfj_jetcharge[ind_rj_rfj_2[1]],key=abs)
        #t_var_comb_jet2_jetChargeMin[0] = min(recojet_rfj_jetcharge[ind_rj_rfj_1[1]],recojet_rfj_jetcharge[ind_rj_rfj_2[1]],key=abs)
        #if recojet_rfj_vector[ind_rj_rfj_1[1]].E()>recojet_rfj_vector[ind_rj_rfj_2[1]].E():
        #    t_var_comb_jet2_E1_jetCharge[0] = recojet_rfj_jetcharge[ind_rj_rfj_1[1]]
        #else:
        #    t_var_comb_jet2_E1_jetCharge[0] = recojet_rfj_jetcharge[ind_rj_rfj_2[1]]

        t_var_comb_jet3_mass[0] = recojet_rfj_vector_combto3[2].M()
        t_var_comb_jet3_theta[0] = recojet_rfj_vector_combto3[2].Theta()
        t_var_comb_jet3_phi[0] = recojet_rfj_vector_combto3[2].Phi()
        t_var_comb_jet3_E[0] = recojet_rfj_vector_combto3[2].E()
        t_var_comb_jet3_BTagMax[0] = max(recojet_rfj_BTag[ind_rj_rfj_1[2]],recojet_rfj_BTag[ind_rj_rfj_2[2]])
        #angles between original jets combined in this one
        t_var_comb_jet3_dalpha[0] = recojet_rfj_vector[ind_rj_rfj_1[2]].Angle(recojet_rfj_vector[ind_rj_rfj_2[2]].Vect())
        t_var_comb_jet3_dphi[0] = recojet_rfj_vector[ind_rj_rfj_1[2]].DeltaPhi(recojet_rfj_vector[ind_rj_rfj_2[2]])
        t_var_comb_jet3_dtheta[0] = abs(recojet_rfj_vector[ind_rj_rfj_1[2]].Theta()-recojet_rfj_vector[ind_rj_rfj_2[2]].Theta())
        #ratio of more energetic input jet for combination to total sum
        t_var_comb_jet3_E1_over_Etot[0] = max(recojet_rfj_vector[ind_rj_rfj_1[2]].E(),recojet_rfj_vector[ind_rj_rfj_2[2]].E())/recojet_rfj_vector_combto3[2].E()
        #t_var_comb_jet3_jetChargeMax[0] = max(recojet_rfj_jetcharge[ind_rj_rfj_1[2]],recojet_rfj_jetcharge[ind_rj_rfj_2[2]],key=abs)
        #t_var_comb_jet3_jetChargeMin[0] = min(recojet_rfj_jetcharge[ind_rj_rfj_1[2]],recojet_rfj_jetcharge[ind_rj_rfj_2[2]],key=abs)
        #if recojet_rfj_vector[ind_rj_rfj_1[2]].E()>recojet_rfj_vector[ind_rj_rfj_2[2]].E():
        #    t_var_comb_jet3_E1_jetCharge[0] = recojet_rfj_jetcharge[ind_rj_rfj_1[2]]
        #else:
        #    t_var_comb_jet3_E1_jetCharge[0] = recojet_rfj_jetcharge[ind_rj_rfj_2[2]]
 
        if len(genjet_vector)>0 and len(genjet_vector_combto3)>0:
            t_var_comb_genjet1_mass[0] = genjet_vector_combto3[0].M()
            t_var_comb_genjet1_theta[0] = genjet_vector_combto3[0].Theta()
            t_var_comb_genjet1_phi[0] = genjet_vector_combto3[0].Phi()
            t_var_comb_genjet1_E[0] = genjet_vector_combto3[0].E()
        #angles between original jets combined in this one
            t_var_comb_genjet1_dalpha[0] = genjet_vector[ind_gj1[0]].Angle(genjet_vector[ind_gj2[0]].Vect())
            t_var_comb_genjet1_dphi[0] = genjet_vector[ind_gj1[0]].DeltaPhi(genjet_vector[ind_gj2[0]])
            t_var_comb_genjet1_dtheta[0] = abs(genjet_vector[ind_gj1[0]].Theta()-genjet_vector[ind_gj2[0]].Theta())
        #ratio of more energetic input jet for combination to total sum
            t_var_comb_genjet1_E1_over_Etot[0] = max(genjet_vector[ind_gj1[0]].E(),genjet_vector[ind_gj2[0]].E())/genjet_vector_combto3[0].E()

        if len(genjet_vector)>1 and len(genjet_vector_combto3)>1:
            t_var_comb_genjet2_mass[0] = genjet_vector_combto3[1].M()
            t_var_comb_genjet2_theta[0] = genjet_vector_combto3[1].Theta()
            t_var_comb_genjet2_phi[0] = genjet_vector_combto3[1].Phi()
            t_var_comb_genjet2_E[0] = genjet_vector_combto3[1].E()
        #angles between original jets combined in this one
            t_var_comb_genjet2_dalpha[0] = genjet_vector[ind_gj1[1]].Angle(genjet_vector[ind_gj2[1]].Vect())
            t_var_comb_genjet2_dphi[0] = genjet_vector[ind_gj1[1]].DeltaPhi(genjet_vector[ind_gj2[1]])
            t_var_comb_genjet2_dtheta[0] = abs(genjet_vector[ind_gj1[1]].Theta()-genjet_vector[ind_gj2[1]].Theta())
        #ratio of more energetic input jet for combination to total sum
            t_var_comb_genjet2_E1_over_Etot[0] = max(genjet_vector[ind_gj1[1]].E(),genjet_vector[ind_gj2[1]].E())/genjet_vector_combto3[1].E()

        if len(genjet_vector)>2 and len(genjet_vector_combto3)>2:
            t_var_comb_genjet3_mass[0] = genjet_vector_combto3[2].M()
            t_var_comb_genjet3_theta[0] = genjet_vector_combto3[2].Theta()
            t_var_comb_genjet3_phi[0] = genjet_vector_combto3[2].Phi()
            t_var_comb_genjet3_E[0] = genjet_vector_combto3[2].E()
        #angles between original jets combined in this one
            t_var_comb_genjet3_dalpha[0] = genjet_vector[ind_gj1[2]].Angle(genjet_vector[ind_gj2[2]].Vect())
            t_var_comb_genjet3_dphi[0] = genjet_vector[ind_gj1[2]].DeltaPhi(genjet_vector[ind_gj2[2]])
            t_var_comb_genjet3_dtheta[0] = abs(genjet_vector[ind_gj1[2]].Theta()-genjet_vector[ind_gj2[2]].Theta())
        #ratio of more energetic input jet for combination to total sum
            t_var_comb_genjet3_E1_over_Etot[0] = max(genjet_vector[ind_gj1[2]].E(),genjet_vector[ind_gj2[2]].E())/genjet_vector_combto3[2].E()

        recojet_rfj_vector_wBTag_E_ordered=recojet_rfj_BTag_Set_sorted[:]
        recojet_rfj_vector_wBTag_E_ordered.sort(key=lambda x: x[0].E(), reverse=True)

        #for rfj in recojet_rfj_vector_wBTag_E_ordered:
        #    print 'rfj with btag ordered after E', rfj[0].E(),rfj[1]
        #for rfj in recojet_rfj_BTag_Set_sorted:
        #    print 'rfj with btag ordered after BTag', rfj[0].E(),rfj[1]     
        #order jets here by energy
        t_var_jet1_E[0] = recojet_rfj_vector_wBTag_E_ordered[0][0].E()
        t_var_jet1_theta[0] = recojet_rfj_vector_wBTag_E_ordered[0][0].Theta()
        t_var_jet1_BTag[0] = recojet_rfj_vector_wBTag_E_ordered[0][1]

        t_var_jet2_E[0] = recojet_rfj_vector_wBTag_E_ordered[1][0].E()
        t_var_jet2_theta[0] = recojet_rfj_vector_wBTag_E_ordered[1][0].Theta()
        t_var_jet2_BTag[0] = recojet_rfj_vector_wBTag_E_ordered[1][1]

        t_var_jet3_E[0] = recojet_rfj_vector_wBTag_E_ordered[2][0].E()
        t_var_jet3_theta[0] = recojet_rfj_vector_wBTag_E_ordered[2][0].Theta()
        t_var_jet3_BTag[0] = recojet_rfj_vector_wBTag_E_ordered[2][1]

        t_var_jet4_E[0] = recojet_rfj_vector_wBTag_E_ordered[3][0].E()
        t_var_jet4_theta[0] = recojet_rfj_vector_wBTag_E_ordered[3][0].Theta()
        t_var_jet4_BTag[0] = recojet_rfj_vector_wBTag_E_ordered[3][1]

        t_var_jet5_E[0] = recojet_rfj_vector_wBTag_E_ordered[4][0].E()
        t_var_jet5_theta[0] = recojet_rfj_vector_wBTag_E_ordered[4][0].Theta()
        t_var_jet5_BTag[0] = recojet_rfj_vector_wBTag_E_ordered[4][1]

        t_var_jet6_E[0] = recojet_rfj_vector_wBTag_E_ordered[5][0].E()
        t_var_jet6_theta[0] = recojet_rfj_vector_wBTag_E_ordered[5][0].Theta()
        t_var_jet6_BTag[0] = recojet_rfj_vector_wBTag_E_ordered[5][1]

        
        t_var_y12[0] = ientry.reco_y12
        t_var_y23[0] = ientry.reco_y23
        t_var_y34[0] = ientry.reco_y34
        t_var_y45[0] = ientry.reco_y45
        t_var_y56[0] = ientry.reco_y56

        t_var_gen_y12[0] = ientry.gen_y12
        t_var_gen_y23[0] = ientry.gen_y23
        t_var_gen_y34[0] = ientry.gen_y34
        t_var_gen_y45[0] = ientry.gen_y45
        t_var_gen_y56[0] = ientry.gen_y56
        
        #order by BTag --> largest BTag is called BTag1
        t_var_BTag1[0] = recojet_rfj_BTag_Set_sorted[0][1]
        t_var_BTag2[0] = recojet_rfj_BTag_Set_sorted[1][1]
        t_var_BTag3[0] = recojet_rfj_BTag_Set_sorted[2][1]
        t_var_BTag4[0] = recojet_rfj_BTag_Set_sorted[3][1]
        
        t_var_BTag_sum_max2[0] = recojet_rfj_BTag_Set_sorted[0][1]+recojet_rfj_BTag_Set_sorted[1][1]
        t_var_BTag_sum_max3[0] = recojet_rfj_BTag_Set_sorted[0][1]+recojet_rfj_BTag_Set_sorted[1][1]+recojet_rfj_BTag_Set_sorted[2][1]
        t_var_BTag_sum_max4[0] = recojet_rfj_BTag_Set_sorted[0][1]+recojet_rfj_BTag_Set_sorted[1][1]+recojet_rfj_BTag_Set_sorted[2][1]+recojet_rfj_BTag_Set_sorted[3][1]
        t_var_BTag_sum_all[0] = recojet_rfj_BTag_Set_sorted[0][1]+recojet_rfj_BTag_Set_sorted[1][1]+recojet_rfj_BTag_Set_sorted[2][1]+recojet_rfj_BTag_Set_sorted[3][1]+recojet_rfj_BTag_Set_sorted[4][1]+recojet_rfj_BTag_Set_sorted[5][1]




        if(recojet_rfj_vector_combto3[0].M()<recojet_rfj_vector_combto3[1].M() or recojet_rfj_vector_combto3[1].M()<recojet_rfj_vector_combto3[2].M() or  recojet_rfj_vector_combto3[0].M()<recojet_rfj_vector_combto3[2].M()):
            print 'mass order screwed up',recojet_rfj_vector_combto3[0].M(),recojet_rfj_vector_combto3[1].M(),recojet_rfj_vector_combto3[1].M()

        totPFO_=TLorentzVector(0,0,0,0)
        totPFO_.SetPxPyPzE(ientry.totPFO_Px,ientry.totPFO_Py,ientry.totPFO_Pz,ientry.totPFO_E)
        tot_trueVis_=TLorentzVector(0,0,0,0)
        tot_trueVis_.SetPxPyPzE(ientry.true_Px,ientry.true_Py,ientry.true_Pz,ientry.true_E)
        tot_trueInv_=TLorentzVector(0,0,0,0)
        tot_trueInv_.SetPxPyPzE(ientry.true_inv_Px,ientry.true_inv_Py,ientry.true_inv_Pz,ientry.true_inv_E)

        t_var_sqrtS[0] =totPFO_.M()
        t_var_sqrtS_gen[0] =tot_trueVis_.M()
        t_var_sqrtS_jets [0] =tot_rfjet_Vector_.M()
        t_var_sqrtS_genjets[0] =tot_genjet_Vector_.M()
        t_var_sqrtS_parton[0]=tempTotEventP4.M()

        t_var_MET[0] =totPFO_.Pt()
        t_var_MET_gen[0] =tot_trueInv_.Pt()
        t_var_MHT[0] =tot_rfjet_Vector_.Pt()
        t_var_MHT_gen[0] =tot_genjet_Vector_.Pt()

        #filled for every entry, after having passed continues for 
        mytree.Fill()
        #for later histo filling
        if len(hist_vec_reco_1D_parton)>0:

        #for vect in range(len(recojet_vector)):
            #    print 'mass vect after ',vect,recojet_vector[vect].M()
            #    a, b = map(list, zip(*y))
            if len(recojet_vector)>2:
                recojet_vector_combto3,ind_rj_1,ind_rj_2=map(list,zip(*orderedThreeVector(recojet_vector)))
                
 


        #for vect in range(len(recojet_vector_combto3)):
        #    print num_entry,'reco mass vect after combination to 3 ',vect,recojet_vector_combto3[vect].M(),recojet_vector_combto3[vect].E(),recojet_vector_combto3[vect].P(),(recojet_vector[ind_rj_1[vect]]+recojet_vector[ind_rj_2[vect]]).M(),(recojet_vector[ind_rj_1[vect]]+recojet_vector[ind_rj_2[vect]]).E()
        

        #print 'index vector reco',len(index_recojet_vector_combto3),index_recojet_vector_combto3
        #print'vector index sum 0',  (recojet_vector[index_recojet_vector_combto3[0]]+recojet_vector[index_recojet_vector_combto3[1]]).M(),(recojet_vector[index_recojet_vector_combto3[0]]+recojet_vector[index_recojet_vector_combto3[1]]).E(),index_group[0]
        #print'vector index sum 1',  (recojet_vector[index_recojet_vector_combto3[2]]+recojet_vector[index_recojet_vector_combto3[3]]).M(),(recojet_vector[index_recojet_vector_combto3[2]]+recojet_vector[index_recojet_vector_combto3[3]]).E(),index_group[1]
        #print'vector index sum 2',  (recojet_vector[index_recojet_vector_combto3[4]]+recojet_vector[index_recojet_vector_combto3[5]]).M(),(recojet_vector[index_recojet_vector_combto3[4]]+recojet_vector[index_recojet_vector_combto3[5]]).E(),index_group[2]

       #checked here that it worked properly, original still not sorted by BTags, set_sorted is indeed properly sorted by BTags
        #for btset in recojet_rfj_BTag_Set:
            #    print 'btagset',btset[0].M(),btset[1]
            
        #for btsetsorted in recojet_rfj_BTag_Set_sorted:
            #    print 'btagset sorted',btsetsorted[0].M(),btsetsorted[1]
            
            case_jet_throwaway=-10
            reco_orderedBTagtuple,case_jet_throwaway=orderedThreeVectorBTagSort(recojet_rfj_BTag_Set_sorted)
            recojet_rfj_vector_combto3_wBTag,ind_rj_rfj_1_btag,ind_rj_rfj_2_btag=map(list,zip(*reco_orderedBTagtuple))
            if len(ind_rj_rfj_1_btag)>len(ind_rj_rfj_2_btag):
                print 'interesting combination, indices are different ', ind_rj_rfj_1_btag,ind_rj_rfj_2_btag, 'in case',case_jet_throwaway
                
                
                for vect in range(len(recojet_rfj_vector_combto3_wBTag)):
                    if ind_rj_rfj_2_btag[vect]!=-1:
                        
                        whateverdosomething=True
                #print num_entry,'reco mass vect after combination to 3 with btag info ',vect,recojet_rfj_vector_combto3_wBTag[vect].M(),recojet_rfj_vector_combto3_wBTag[vect].E(),recojet_rfj_vector_combto3_wBTag[vect].P(),ind_rj_rfj_1_btag[vect],ind_rj_rfj_2_btag[vect],recojet_rfj_BTag_Set_sorted[ind_rj_rfj_1_btag[vect]][1],recojet_rfj_BTag_Set_sorted[ind_rj_rfj_2_btag[vect]][1]
                    else:
                        whateverdosomething=False
                #print num_entry,'reco mass vect after combination to 3 with btag info ',vect,recojet_rfj_vector_combto3_wBTag[vect].M(),recojet_rfj_vector_combto3_wBTag[vect].E(),recojet_rfj_vector_combto3_wBTag[vect].P(),ind_rj_rfj_1_btag[vect],ind_rj_rfj_2_btag[vect],recojet_rfj_BTag_Set_sorted[ind_rj_rfj_1_btag[vect]][1]
        #for vect in recojet_rfj_BTag_Set_sorted:
        #    print 'mass vect after BTag sort',vect[0].M(),vect[1]
        #for vect in range(len(recojet_vector)):
        #    print 'mass vect after ',vect,recojet_vector[vect].M()
                        


            hist_vec_reco_1D_parton[113].Fill(temp_gj_sum.M(),weight)
            hist_vec_reco_1D_parton[115].Fill((temp_true_inv+temp_gj_sum).M(),weight)
            hist_vec_reco_1D_parton[117].Fill(tempTotEventP4.M(),weight)
            for ind in range(len(recojet_rfj_E)):
                temp_rfj_=TLorentzVector(0,0,0,0)
                temp_rfj_.SetPxPyPzE(recojet_rfj_Px[ind],recojet_rfj_Py[ind],recojet_rfj_Pz[ind],recojet_rfj_E[ind])
            #check what's going on with high mass jets
                if temp_rfj_.M()>150:
                #print 'very high value here',temp_rfj_.M()
                    if recojet_rfj_BTag[ind]<0.5:
                        hist_vec_reco_1D_parton[190].Fill(degrees(min(temp_rfj_.Angle(tempH1P4.Vect()),temp_rfj_.Angle(tempH2P4.Vect()))),weight)
                        hist_vec_reco_1D_parton[191].Fill(degrees(temp_rfj_.Angle(tempZ_first.Vect())),weight)
                        hist_vec_reco_1D_parton[192].Fill(degrees(min(temp_rfj_.Angle(tempH1_q2.Vect()),temp_rfj_.Angle(tempH1_q1.Vect()),temp_rfj_.Angle(tempH2_q2.Vect()),temp_rfj_.Angle(tempH2_q1.Vect()))),weight)
                        hist_vec_reco_1D_parton[193].Fill(degrees(min(temp_rfj_.Angle(tempZ_q_pos.Vect()),temp_rfj_.Angle(tempZ_q_neg.Vect()))),weight)
                        if len(genjet_vector)>5:
                            hist_vec_reco_1D_parton[198].Fill(degrees(min(temp_rfj_.Angle(genjet_vector[0].Vect()),temp_rfj_.Angle(genjet_vector[1].Vect()),temp_rfj_.Angle(genjet_vector[2].Vect()),temp_rfj_.Angle(genjet_vector[3].Vect()),temp_rfj_.Angle(genjet_vector[4].Vect()),temp_rfj_.Angle(genjet_vector[5].Vect()))),weight)
                    else:
                        hist_vec_reco_1D_parton[194].Fill(degrees(min(temp_rfj_.Angle(tempH1P4.Vect()),temp_rfj_.Angle(tempH2P4.Vect()))),weight)
                        hist_vec_reco_1D_parton[195].Fill(degrees(temp_rfj_.Angle(tempZ_first.Vect())),weight)
                        hist_vec_reco_1D_parton[196].Fill(degrees(min(temp_rfj_.Angle(tempH1_q2.Vect()),temp_rfj_.Angle(tempH1_q1.Vect()),temp_rfj_.Angle(tempH2_q2.Vect()),temp_rfj_.Angle(tempH2_q1.Vect()))),weight)
                        hist_vec_reco_1D_parton[197].Fill(degrees(min(temp_rfj_.Angle(tempZ_q_pos.Vect()),temp_rfj_.Angle(tempZ_q_neg.Vect()))),weight)
                        if len(genjet_vector)>5:
                            hist_vec_reco_1D_parton[199].Fill(degrees(min(temp_rfj_.Angle(genjet_vector[0].Vect()),temp_rfj_.Angle(genjet_vector[1].Vect()),temp_rfj_.Angle(genjet_vector[2].Vect()),temp_rfj_.Angle(genjet_vector[3].Vect()),temp_rfj_.Angle(genjet_vector[4].Vect()),temp_rfj_.Angle(genjet_vector[5].Vect()))),weight)


            if(ind_rj_rfj_2_btag[0]!=-1):
                hist_vec_reco_1D_parton[221].Fill(max(recojet_rfj_BTag_Set_sorted[ind_rj_rfj_1_btag[0]][1],recojet_rfj_BTag_Set_sorted[ind_rj_rfj_2_btag[0]][1]),weight)
            else:
                hist_vec_reco_1D_parton[221].Fill(recojet_rfj_BTag_Set_sorted[ind_rj_rfj_1_btag[0]][1],weight)
            if(ind_rj_rfj_2_btag[1]!=-1):
                hist_vec_reco_1D_parton[222].Fill(max(recojet_rfj_BTag_Set_sorted[ind_rj_rfj_1_btag[1]][1],recojet_rfj_BTag_Set_sorted[ind_rj_rfj_2_btag[1]][1]),weight)
            else:
                hist_vec_reco_1D_parton[222].Fill(recojet_rfj_BTag_Set_sorted[ind_rj_rfj_1_btag[1]][1],weight)
            if(ind_rj_rfj_2_btag[2]!=-1):
                hist_vec_reco_1D_parton[223].Fill(max(recojet_rfj_BTag_Set_sorted[ind_rj_rfj_1_btag[2]][1],recojet_rfj_BTag_Set_sorted[ind_rj_rfj_2_btag[2]][1]),weight)
            else:
                hist_vec_reco_1D_parton[223].Fill(recojet_rfj_BTag_Set_sorted[ind_rj_rfj_1_btag[2]][1],weight)

            hist_vec_reco_1D_parton[206].Fill(recojet_rfj_vector_combto3_wBTag[0].M(),weight)
            hist_vec_reco_1D_parton[207].Fill(recojet_rfj_vector_combto3_wBTag[1].M(),weight)
            hist_vec_reco_1D_parton[208].Fill(recojet_rfj_vector_combto3_wBTag[2].M(),weight)
            
            hist_vec_reco_1D_parton[209].Fill(recojet_rfj_vector_combto3_wBTag[0].M(),weight)
            hist_vec_reco_1D_parton[210].Fill(recojet_rfj_vector_combto3_wBTag[1].M(),weight)
            hist_vec_reco_1D_parton[211].Fill(recojet_rfj_vector_combto3_wBTag[2].M(),weight)


        #first comb_wBTag rfj jet and H1
            if(recojet_rfj_vector_combto3_wBTag[0].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3_wBTag[0].Angle(tempH2P4.Vect()) and recojet_rfj_vector_combto3_wBTag[0].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3_wBTag[0].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[230].Fill(recojet_rfj_vector_combto3_wBTag[0].P()/tempH1P4.P(),weight)
       #first comb_wBTag rfj jet and H2
            elif(recojet_rfj_vector_combto3_wBTag[0].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3_wBTag[0].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3_wBTag[0].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3_wBTag[0].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[230].Fill(recojet_rfj_vector_combto3_wBTag[0].P()/tempH2P4.P(),weight)
      #first comb_wBTag rfj jet and Z
            elif(recojet_rfj_vector_combto3_wBTag[0].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3_wBTag[0].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3_wBTag[0].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3_wBTag[0].Angle(tempH2P4.Vect())):
                hist_vec_reco_1D_parton[230].Fill(recojet_rfj_vector_combto3_wBTag[0].P()/tempZP4.P(),weight)
            else:
                print 'btag rfj rj1 should have done all combinations'

        #second comb_wBTag rfj jet and H1
            if(recojet_rfj_vector_combto3_wBTag[1].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3_wBTag[1].Angle(tempH2P4.Vect()) and recojet_rfj_vector_combto3_wBTag[1].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3_wBTag[1].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[231].Fill(recojet_rfj_vector_combto3_wBTag[1].P()/tempH1P4.P(),weight)
       #second comb_wBTag rfj jet and H2
            elif(recojet_rfj_vector_combto3_wBTag[1].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3_wBTag[1].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3_wBTag[1].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3_wBTag[1].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[231].Fill(recojet_rfj_vector_combto3_wBTag[1].P()/tempH2P4.P(),weight)
      #second comb_wBTag rfj jet and Z
            elif(recojet_rfj_vector_combto3_wBTag[1].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3_wBTag[1].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3_wBTag[1].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3_wBTag[1].Angle(tempH2P4.Vect())):
                hist_vec_reco_1D_parton[231].Fill(recojet_rfj_vector_combto3_wBTag[2].P()/tempZP4.P(),weight)
            else:
                print 'btag rfj rj2 should have done all combinations'

        #third comb_wBTag rfj jet and H1
            if(recojet_rfj_vector_combto3_wBTag[2].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3_wBTag[2].Angle(tempH2P4.Vect()) and recojet_rfj_vector_combto3_wBTag[2].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3_wBTag[2].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[232].Fill(recojet_rfj_vector_combto3_wBTag[2].P()/tempH1P4.P(),weight)
        #third comb_wBTag rfj jet and H2
            elif(recojet_rfj_vector_combto3_wBTag[2].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3_wBTag[2].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3_wBTag[2].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3_wBTag[2].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[232].Fill(recojet_rfj_vector_combto3_wBTag[2].P()/tempH2P4.P(),weight)
        #third comb_wBTag rfj jet and Z
            elif(recojet_rfj_vector_combto3_wBTag[2].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3_wBTag[2].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3_wBTag[2].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3_wBTag[2].Angle(tempH2P4.Vect())):
                hist_vec_reco_1D_parton[232].Fill(recojet_rfj_vector_combto3_wBTag[2].P()/tempZP4.P(),weight)
            else:
                print 'btag rfj rj3 should have done all combinations'



            for rj_wbtag in range(len(recojet_rfj_vector_combto3_wBTag)):
                min_angle = 200.
                if len(genjet_vector)>2: 
                    for gj in genjet_vector_combto3:
                        if( degrees(recojet_rfj_vector_combto3_wBTag[rj_wbtag].Angle(gj.Vect()))<min_angle) :
                            min_angle=degrees(recojet_rfj_vector_combto3_wBTag[rj_wbtag].Angle(gj.Vect()))
                    if rj_wbtag==0:
                        hist_vec_reco_1D_parton[212].Fill(min_angle,weight)
                    elif rj_wbtag==1:
                        hist_vec_reco_1D_parton[213].Fill(min_angle,weight)
                    elif rj_wbtag==2:
                        hist_vec_reco_1D_parton[214].Fill(min_angle,weight)

            hist_vec_reco_1D_parton[215].Fill(degrees(min(recojet_rfj_vector_combto3_wBTag[0].Angle(tempH1P4.Vect()),recojet_rfj_vector_combto3_wBTag[0].Angle(tempH2P4.Vect()))),weight)
            hist_vec_reco_1D_parton[216].Fill(degrees(recojet_rfj_vector_combto3_wBTag[0].Angle(tempZP4.Vect())),weight)
            hist_vec_reco_1D_parton[217].Fill(degrees(min(recojet_rfj_vector_combto3_wBTag[1].Angle(tempH1P4.Vect()),recojet_rfj_vector_combto3_wBTag[1].Angle(tempH2P4.Vect()))),weight)
            hist_vec_reco_1D_parton[218].Fill(degrees(recojet_rfj_vector_combto3_wBTag[1].Angle(tempZP4.Vect())),weight)
            hist_vec_reco_1D_parton[219].Fill(degrees(min(recojet_rfj_vector_combto3_wBTag[2].Angle(tempH1P4.Vect()),recojet_rfj_vector_combto3_wBTag[2].Angle(tempH2P4.Vect()))),weight)
            hist_vec_reco_1D_parton[220].Fill(degrees(recojet_rfj_vector_combto3_wBTag[2].Angle(tempZP4.Vect())),weight)


       #now histos of H1 and combination of rfj vectors
            if(tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect()) and tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())):
                hist_vec_reco_1D_parton[239].Fill(recojet_rfj_vector_combto3_wBTag[0].P()/tempH1P4.P(),weight)
                hist_vec_reco_1D_parton[248].Fill(degrees(tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect())),weight)
            elif(tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect()) and tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())):
                hist_vec_reco_1D_parton[239].Fill(recojet_rfj_vector_combto3_wBTag[1].P()/tempH1P4.P(),weight)
                hist_vec_reco_1D_parton[248].Fill(degrees(tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())),weight)
            elif(tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect()) and tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())):
                hist_vec_reco_1D_parton[239].Fill(recojet_rfj_vector_combto3_wBTag[2].P()/tempH1P4.P(),weight)
                hist_vec_reco_1D_parton[248].Fill(degrees(tempH1P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())),weight)
            else:
                print 'sth wrong in H1 to rfj rj com wBTag lineup'

       #now histos of H2 and combination of rfj vectors
            if(tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect()) and tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())):
                hist_vec_reco_1D_parton[240].Fill(recojet_rfj_vector_combto3_wBTag[0].P()/tempH2P4.P(),weight)
                hist_vec_reco_1D_parton[249].Fill(degrees(tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect())),weight)
            elif(tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect()) and tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())):
                hist_vec_reco_1D_parton[240].Fill(recojet_rfj_vector_combto3_wBTag[1].P()/tempH2P4.P(),weight)
                hist_vec_reco_1D_parton[249].Fill(degrees(tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())),weight)
            elif(tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect()) and tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())):
                hist_vec_reco_1D_parton[240].Fill(recojet_rfj_vector_combto3_wBTag[2].P()/tempH2P4.P(),weight)
                hist_vec_reco_1D_parton[249].Fill(degrees(tempH2P4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())),weight)
            else:
                print 'sth wrong in H2 to rfj rj com wBTag lineup'

      #now histos of Z and combination of rfj vectors
            if(tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect()) and tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())):
                hist_vec_reco_1D_parton[241].Fill(recojet_rfj_vector_combto3_wBTag[0].P()/tempZP4.P(),weight)
                hist_vec_reco_1D_parton[250].Fill(degrees(tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect())),weight)
            elif(tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect()) and tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())):
                hist_vec_reco_1D_parton[241].Fill(recojet_rfj_vector_combto3_wBTag[1].P()/tempZP4.P(),weight)
                hist_vec_reco_1D_parton[250].Fill(degrees(tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())),weight)
            elif(tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[0].Vect()) and tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[1].Vect())):
                hist_vec_reco_1D_parton[241].Fill(recojet_rfj_vector_combto3_wBTag[2].P()/tempZP4.P(),weight)
                hist_vec_reco_1D_parton[250].Fill(degrees(tempZP4.Angle(recojet_rfj_vector_combto3_wBTag[2].Vect())),weight)
            else:
                print 'sth wrong in Z to rfj rj com wBTag lineup'

            if len(recojet_vector)>2:
                if(tempH1P4.E()>tempH2P4.E()):
                    hist_vec_reco_1D_parton[77].Fill(degrees(recojet_vector[0].Angle(tempH1P4.Vect())),weight)
                    hist_vec_reco_1D_parton[78].Fill(degrees(recojet_vector[0].Angle(tempH2P4.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[80].Fill(degrees(recojet_vector[1].Angle(tempH1P4.Vect())),weight)
                    hist_vec_reco_1D_parton[81].Fill(degrees(recojet_vector[1].Angle(tempH2P4.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[83].Fill(degrees(recojet_vector[2].Angle(tempH1P4.Vect())),weight)
                    hist_vec_reco_1D_parton[84].Fill(degrees(recojet_vector[2].Angle(tempH2P4.Vect())),weight)
                else:
                    hist_vec_reco_1D_parton[77].Fill(degrees(recojet_vector[0].Angle(tempH2P4.Vect())),weight)
                    hist_vec_reco_1D_parton[78].Fill(degrees(recojet_vector[0].Angle(tempH1P4.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[80].Fill(degrees(recojet_vector[1].Angle(tempH2P4.Vect())),weight)
                    hist_vec_reco_1D_parton[81].Fill(degrees(recojet_vector[1].Angle(tempH1P4.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[83].Fill(degrees(recojet_vector[2].Angle(tempH2P4.Vect())),weight)
                    hist_vec_reco_1D_parton[84].Fill(degrees(recojet_vector[2].Angle(tempH1P4.Vect())),weight)

                hist_vec_reco_1D_parton[79].Fill(degrees(recojet_vector[1].Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[82].Fill(degrees(recojet_vector[1].Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[85].Fill(degrees(recojet_vector[2].Angle(tempZP4.Vect())),weight)
                
                hist_vec_reco_1D_parton[86].Fill(recojet_vector[0].M(),weight)
                hist_vec_reco_1D_parton[87].Fill(recojet_vector[1].M(),weight)
                hist_vec_reco_1D_parton[88].Fill(recojet_vector[2].M(),weight)
                hist_vec_reco_1D_parton[133].Fill(recojet_vector[0].E(),weight)
                hist_vec_reco_1D_parton[134].Fill(recojet_vector[1].E(),weight)
                hist_vec_reco_1D_parton[135].Fill(recojet_vector[2].E(),weight)
            if len(recojet_vector)>3:
                hist_vec_reco_1D_parton[127].Fill(degrees(min(recojet_vector[3].Angle(tempH1P4.Vect()),recojet_vector[3].Angle(tempH2P4.Vect()))),weight)
                hist_vec_reco_1D_parton[128].Fill(degrees(recojet_vector[3].Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[124].Fill(recojet_vector[3].M(),weight)
                hist_vec_reco_1D_parton[136].Fill(recojet_vector[2].E(),weight)
                if len(recojet_vector)>4:
                    hist_vec_reco_1D_parton[139].Fill(degrees(min(recojet_vector[4].Angle(tempH1P4.Vect()),recojet_vector[4].Angle(tempH2P4.Vect()))),weight)
                    hist_vec_reco_1D_parton[140].Fill(degrees(recojet_vector[4].Angle(tempZP4.Vect())),weight)
                    if len(recojet_vector)>5:
                        hist_vec_reco_1D_parton[143].Fill(degrees(min(recojet_vector[5].Angle(tempH1P4.Vect()),recojet_vector[5].Angle(tempH2P4.Vect()))),weight)
                        hist_vec_reco_1D_parton[144].Fill(degrees(recojet_vector[5].Angle(tempZP4.Vect())),weight)


            for rj in range(len(recojet_vector)):
                min_angle = 200.
                if len(genjet_vector)>2: 
                   for gj in genjet_vector:
                       if( degrees(recojet_vector[rj].Angle(gj.Vect()))<min_angle) :
                           min_angle=degrees(recojet_vector[rj].Angle(gj.Vect()))
                   if rj==0:
                       hist_vec_reco_1D_parton[145].Fill(min_angle,weight)
                   elif rj==1:
                       hist_vec_reco_1D_parton[146].Fill(min_angle,weight)
                   elif rj==2:
                       hist_vec_reco_1D_parton[147].Fill(min_angle,weight)
                   elif rj==3:
                       hist_vec_reco_1D_parton[148].Fill(min_angle,weight)
                   elif rj==4:
                       hist_vec_reco_1D_parton[149].Fill(min_angle,weight)
                   elif rj==5:
                       hist_vec_reco_1D_parton[150].Fill(min_angle,weight)

            hist_vec_reco_1D_parton[114].Fill(temp_gj_sum.M(),weight)
            hist_vec_reco_1D_parton[116].Fill((temp_true_inv+temp_gj_sum).M(),weight)
            if len(recojet_vector)>2:
                hist_vec_reco_1D_parton[154].Fill(recojet_vector_combto3[0].M(),weight)
                hist_vec_reco_1D_parton[155].Fill(recojet_vector_combto3[1].M(),weight)
                hist_vec_reco_1D_parton[156].Fill(recojet_vector_combto3[2].M(),weight)
                
                hist_vec_reco_1D_parton[200].Fill(recojet_vector_combto3[0].M(),weight)
                hist_vec_reco_1D_parton[201].Fill(recojet_vector_combto3[1].M(),weight)
                hist_vec_reco_1D_parton[202].Fill(recojet_vector_combto3[2].M(),weight)

        #first comb jet and H1
                if(recojet_vector_combto3[0].Angle(tempH1P4.Vect())<recojet_vector_combto3[0].Angle(tempH2P4.Vect()) and recojet_vector_combto3[0].Angle(tempH1P4.Vect())<recojet_vector_combto3[0].Angle(tempZP4.Vect())):
                    hist_vec_reco_1D_parton[224].Fill(recojet_vector_combto3[0].P()/tempH1P4.P(),weight)
       #first comb jet and H2
                elif(recojet_vector_combto3[0].Angle(tempH2P4.Vect())<recojet_vector_combto3[0].Angle(tempH1P4.Vect()) and recojet_vector_combto3[0].Angle(tempH2P4.Vect())<recojet_vector_combto3[0].Angle(tempZP4.Vect())):
                    hist_vec_reco_1D_parton[224].Fill(recojet_vector_combto3[0].P()/tempH2P4.P(),weight)
      #first comb jet and Z
                elif(recojet_vector_combto3[0].Angle(tempZP4.Vect())<recojet_vector_combto3[0].Angle(tempH1P4.Vect()) and recojet_vector_combto3[0].Angle(tempZP4.Vect())<recojet_vector_combto3[0].Angle(tempH2P4.Vect())):
                    hist_vec_reco_1D_parton[224].Fill(recojet_vector_combto3[0].P()/tempZP4.P(),weight)
                else:
                    print 'rj1 should have done all combinations'

        #second comb jet and H1
                if(recojet_vector_combto3[1].Angle(tempH1P4.Vect())<recojet_vector_combto3[1].Angle(tempH2P4.Vect()) and recojet_vector_combto3[1].Angle(tempH1P4.Vect())<recojet_vector_combto3[1].Angle(tempZP4.Vect())):
                    hist_vec_reco_1D_parton[225].Fill(recojet_vector_combto3[1].P()/tempH1P4.P(),weight)
       #second comb jet and H2
                elif(recojet_vector_combto3[1].Angle(tempH2P4.Vect())<recojet_vector_combto3[1].Angle(tempH1P4.Vect()) and recojet_vector_combto3[1].Angle(tempH2P4.Vect())<recojet_vector_combto3[1].Angle(tempZP4.Vect())):
                    hist_vec_reco_1D_parton[225].Fill(recojet_vector_combto3[1].P()/tempH2P4.P(),weight)
      #second comb jet and Z
                elif(recojet_vector_combto3[1].Angle(tempZP4.Vect())<recojet_vector_combto3[1].Angle(tempH1P4.Vect()) and recojet_vector_combto3[1].Angle(tempZP4.Vect())<recojet_vector_combto3[1].Angle(tempH2P4.Vect())):
                    hist_vec_reco_1D_parton[225].Fill(recojet_vector_combto3[2].P()/tempZP4.P(),weight)
                else:
                    print 'rj2 should have done all combinations'

        #third comb jet and H1
                if(recojet_vector_combto3[2].Angle(tempH1P4.Vect())<recojet_vector_combto3[2].Angle(tempH2P4.Vect()) and recojet_vector_combto3[2].Angle(tempH1P4.Vect())<recojet_vector_combto3[2].Angle(tempZP4.Vect())):
                    hist_vec_reco_1D_parton[226].Fill(recojet_vector_combto3[2].P()/tempH1P4.P(),weight)
        #third comb jet and H2
                elif(recojet_vector_combto3[2].Angle(tempH2P4.Vect())<recojet_vector_combto3[2].Angle(tempH1P4.Vect()) and recojet_vector_combto3[2].Angle(tempH2P4.Vect())<recojet_vector_combto3[2].Angle(tempZP4.Vect())):
                    hist_vec_reco_1D_parton[226].Fill(recojet_vector_combto3[2].P()/tempH2P4.P(),weight)
        #third comb jet and Z
                elif(recojet_vector_combto3[2].Angle(tempZP4.Vect())<recojet_vector_combto3[2].Angle(tempH1P4.Vect()) and recojet_vector_combto3[2].Angle(tempZP4.Vect())<recojet_vector_combto3[2].Angle(tempH2P4.Vect())):
                    hist_vec_reco_1D_parton[226].Fill(recojet_vector_combto3[2].P()/tempZP4.P(),weight)
                else:
                    print 'rj3 should have done all combinations'
                
        #now histos of H1 and combination of vectors
                if(tempH1P4.Angle(recojet_vector_combto3[0].Vect())<tempH1P4.Angle(recojet_vector_combto3[1].Vect()) and tempH1P4.Angle(recojet_vector_combto3[0].Vect())<tempH1P4.Angle(recojet_vector_combto3[2].Vect())):
                    hist_vec_reco_1D_parton[233].Fill(recojet_vector_combto3[0].P()/tempH1P4.P(),weight)
                    hist_vec_reco_1D_parton[242].Fill(degrees(tempH1P4.Angle(recojet_vector_combto3[0].Vect())),weight)
                elif(tempH1P4.Angle(recojet_vector_combto3[1].Vect())<tempH1P4.Angle(recojet_vector_combto3[0].Vect()) and tempH1P4.Angle(recojet_vector_combto3[1].Vect())<tempH1P4.Angle(recojet_vector_combto3[2].Vect())):
                    hist_vec_reco_1D_parton[233].Fill(recojet_vector_combto3[1].P()/tempH1P4.P(),weight)
                    hist_vec_reco_1D_parton[242].Fill(degrees(tempH1P4.Angle(recojet_vector_combto3[1].Vect())),weight)
                elif(tempH1P4.Angle(recojet_vector_combto3[2].Vect())<tempH1P4.Angle(recojet_vector_combto3[0].Vect()) and tempH1P4.Angle(recojet_vector_combto3[2].Vect())<tempH1P4.Angle(recojet_vector_combto3[1].Vect())):
                    hist_vec_reco_1D_parton[233].Fill(recojet_vector_combto3[2].P()/tempH1P4.P(),weight)
                    hist_vec_reco_1D_parton[242].Fill(degrees(tempH1P4.Angle(recojet_vector_combto3[2].Vect())),weight)
                else:
                    print 'sth wrong in H1 to rj comb lineup'

       #now histos of H2 and combination of vectors
                if(tempH2P4.Angle(recojet_vector_combto3[0].Vect())<tempH2P4.Angle(recojet_vector_combto3[1].Vect()) and tempH2P4.Angle(recojet_vector_combto3[0].Vect())<tempH2P4.Angle(recojet_vector_combto3[2].Vect())):
                    hist_vec_reco_1D_parton[234].Fill(recojet_vector_combto3[0].P()/tempH2P4.P(),weight)
                    hist_vec_reco_1D_parton[243].Fill(degrees(tempH2P4.Angle(recojet_vector_combto3[0].Vect())),weight)
                elif(tempH2P4.Angle(recojet_vector_combto3[1].Vect())<tempH2P4.Angle(recojet_vector_combto3[0].Vect()) and tempH2P4.Angle(recojet_vector_combto3[1].Vect())<tempH2P4.Angle(recojet_vector_combto3[2].Vect())):
                    hist_vec_reco_1D_parton[234].Fill(recojet_vector_combto3[1].P()/tempH2P4.P(),weight)
                    hist_vec_reco_1D_parton[243].Fill(degrees(tempH2P4.Angle(recojet_vector_combto3[1].Vect())),weight)
                elif(tempH2P4.Angle(recojet_vector_combto3[2].Vect())<tempH2P4.Angle(recojet_vector_combto3[0].Vect()) and tempH2P4.Angle(recojet_vector_combto3[2].Vect())<tempH2P4.Angle(recojet_vector_combto3[1].Vect())):
                    hist_vec_reco_1D_parton[234].Fill(recojet_vector_combto3[2].P()/tempH2P4.P(),weight)
                    hist_vec_reco_1D_parton[243].Fill(degrees(tempH2P4.Angle(recojet_vector_combto3[2].Vect())),weight)
                else:
                    print 'sth wrong in H2 to rj comb lineup'
                
      #now histos of Z and combination of vectors
                if(tempZP4.Angle(recojet_vector_combto3[0].Vect())<tempZP4.Angle(recojet_vector_combto3[1].Vect()) and tempZP4.Angle(recojet_vector_combto3[0].Vect())<tempZP4.Angle(recojet_vector_combto3[2].Vect())):
                    hist_vec_reco_1D_parton[235].Fill(recojet_vector_combto3[0].P()/tempZP4.P(),weight)
                    hist_vec_reco_1D_parton[244].Fill(degrees(tempZP4.Angle(recojet_vector_combto3[0].Vect())),weight)
                elif(tempZP4.Angle(recojet_vector_combto3[1].Vect())<tempZP4.Angle(recojet_vector_combto3[0].Vect()) and tempZP4.Angle(recojet_vector_combto3[1].Vect())<tempZP4.Angle(recojet_vector_combto3[2].Vect())):
                    hist_vec_reco_1D_parton[235].Fill(recojet_vector_combto3[1].P()/tempZP4.P(),weight)
                    hist_vec_reco_1D_parton[244].Fill(degrees(tempZP4.Angle(recojet_vector_combto3[1].Vect())),weight)
                elif(tempZP4.Angle(recojet_vector_combto3[2].Vect())<tempZP4.Angle(recojet_vector_combto3[0].Vect()) and tempZP4.Angle(recojet_vector_combto3[2].Vect())<tempZP4.Angle(recojet_vector_combto3[1].Vect())):
                    hist_vec_reco_1D_parton[235].Fill(recojet_vector_combto3[2].P()/tempZP4.P(),weight)
                    hist_vec_reco_1D_parton[244].Fill(degrees(tempZP4.Angle(recojet_vector_combto3[2].Vect())),weight)
                else:
                    print 'sth wrong in Z to rj comb lineup'

                for rj in range(len(recojet_vector_combto3)):
                    min_angle = 200.
                    if len(genjet_vector)>2: 
                        for gj in genjet_vector_combto3:
                            if( degrees(recojet_vector_combto3[rj].Angle(gj.Vect()))<min_angle) :
                                min_angle=degrees(recojet_vector_combto3[rj].Angle(gj.Vect()))
                        if rj==0:
                            hist_vec_reco_1D_parton[157].Fill(min_angle,weight)
                        elif rj==1:
                            hist_vec_reco_1D_parton[158].Fill(min_angle,weight)
                        elif rj==2:
                            hist_vec_reco_1D_parton[159].Fill(min_angle,weight)
                    
                hist_vec_reco_1D_parton[160].Fill(degrees(min(recojet_vector_combto3[0].Angle(tempH1P4.Vect()),recojet_vector_combto3[0].Angle(tempH2P4.Vect()))),weight)
                hist_vec_reco_1D_parton[161].Fill(degrees(recojet_vector_combto3[0].Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[162].Fill(degrees(min(recojet_vector_combto3[1].Angle(tempH1P4.Vect()),recojet_vector_combto3[1].Angle(tempH2P4.Vect()))),weight)
                hist_vec_reco_1D_parton[163].Fill(degrees(recojet_vector_combto3[1].Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[164].Fill(degrees(min(recojet_vector_combto3[2].Angle(tempH1P4.Vect()),recojet_vector_combto3[2].Angle(tempH2P4.Vect()))),weight)
                hist_vec_reco_1D_parton[165].Fill(degrees(recojet_vector_combto3[2].Angle(tempZP4.Vect())),weight)
            
 
            hist_vec_reco_1D_parton[178].Fill(recojet_rfj_BTag[ind_rj_rfj_1[0]]+recojet_rfj_BTag[ind_rj_rfj_2[0]],weight)
            hist_vec_reco_1D_parton[179].Fill(recojet_rfj_BTag[ind_rj_rfj_1[1]]+recojet_rfj_BTag[ind_rj_rfj_2[1]],weight)
            hist_vec_reco_1D_parton[180].Fill(recojet_rfj_BTag[ind_rj_rfj_1[2]]+recojet_rfj_BTag[ind_rj_rfj_2[2]],weight)
            hist_vec_reco_1D_parton[181].Fill(max(recojet_rfj_BTag[ind_rj_rfj_1[0]],recojet_rfj_BTag[ind_rj_rfj_2[0]]),weight)
            hist_vec_reco_1D_parton[182].Fill(max(recojet_rfj_BTag[ind_rj_rfj_1[1]],recojet_rfj_BTag[ind_rj_rfj_2[1]]),weight)
            hist_vec_reco_1D_parton[183].Fill(max(recojet_rfj_BTag[ind_rj_rfj_1[2]],recojet_rfj_BTag[ind_rj_rfj_2[2]]),weight)
            
            
            hist_vec_reco_1D_parton[166].Fill(recojet_rfj_vector_combto3[0].M(),weight)
            hist_vec_reco_1D_parton[167].Fill(recojet_rfj_vector_combto3[1].M(),weight)
            hist_vec_reco_1D_parton[168].Fill(recojet_rfj_vector_combto3[2].M(),weight)
            
            hist_vec_reco_1D_parton[203].Fill(recojet_rfj_vector_combto3[0].M(),weight)
            hist_vec_reco_1D_parton[204].Fill(recojet_rfj_vector_combto3[1].M(),weight)
            hist_vec_reco_1D_parton[205].Fill(recojet_rfj_vector_combto3[2].M(),weight)

            for rj in range(len(recojet_rfj_vector_combto3)):
                min_angle = 200.
                if len(genjet_vector)>2: 
                    for gj in genjet_vector_combto3:
                        if( degrees(recojet_rfj_vector_combto3[rj].Angle(gj.Vect()))<min_angle) :
                            min_angle=degrees(recojet_rfj_vector_combto3[rj].Angle(gj.Vect()))
                    if rj==0:
                        hist_vec_reco_1D_parton[169].Fill(min_angle,weight)
                    elif rj==1:
                        hist_vec_reco_1D_parton[170].Fill(min_angle,weight)
                    elif rj==2:
                        hist_vec_reco_1D_parton[171].Fill(min_angle,weight)

            hist_vec_reco_1D_parton[172].Fill(degrees(min(recojet_rfj_vector_combto3[0].Angle(tempH1P4.Vect()),recojet_rfj_vector_combto3[0].Angle(tempH2P4.Vect()))),weight)
            hist_vec_reco_1D_parton[173].Fill(degrees(recojet_rfj_vector_combto3[0].Angle(tempZP4.Vect())),weight)
            hist_vec_reco_1D_parton[174].Fill(degrees(min(recojet_rfj_vector_combto3[1].Angle(tempH1P4.Vect()),recojet_rfj_vector_combto3[1].Angle(tempH2P4.Vect()))),weight)
            hist_vec_reco_1D_parton[175].Fill(degrees(recojet_rfj_vector_combto3[1].Angle(tempZP4.Vect())),weight)
            hist_vec_reco_1D_parton[176].Fill(degrees(min(recojet_rfj_vector_combto3[2].Angle(tempH1P4.Vect()),recojet_rfj_vector_combto3[2].Angle(tempH2P4.Vect()))),weight)
            hist_vec_reco_1D_parton[177].Fill(degrees(recojet_rfj_vector_combto3[2].Angle(tempZP4.Vect())),weight)

        #first comb rfj jet and H1
            if(recojet_rfj_vector_combto3[0].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3[0].Angle(tempH2P4.Vect()) and recojet_rfj_vector_combto3[0].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3[0].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[227].Fill(recojet_rfj_vector_combto3[0].P()/tempH1P4.P(),weight)
       #first comb rfj jet and H2
            elif(recojet_rfj_vector_combto3[0].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3[0].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3[0].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3[0].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[227].Fill(recojet_rfj_vector_combto3[0].P()/tempH2P4.P(),weight)
      #first comb rfj jet and Z
            elif(recojet_rfj_vector_combto3[0].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3[0].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3[0].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3[0].Angle(tempH2P4.Vect())):
                hist_vec_reco_1D_parton[227].Fill(recojet_rfj_vector_combto3[0].P()/tempZP4.P(),weight)
            else:
                print 'rfj rj1 should have done all combinations'

        #second comb rfj jet and H1
            if(recojet_rfj_vector_combto3[1].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3[1].Angle(tempH2P4.Vect()) and recojet_rfj_vector_combto3[1].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3[1].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[228].Fill(recojet_rfj_vector_combto3[1].P()/tempH1P4.P(),weight)
       #second comb rfj jet and H2
            elif(recojet_rfj_vector_combto3[1].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3[1].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3[1].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3[1].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[228].Fill(recojet_rfj_vector_combto3[1].P()/tempH2P4.P(),weight)
      #second comb rfj jet and Z
            elif(recojet_rfj_vector_combto3[1].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3[1].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3[1].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3[1].Angle(tempH2P4.Vect())):
                hist_vec_reco_1D_parton[228].Fill(recojet_rfj_vector_combto3[2].P()/tempZP4.P(),weight)
            else:
                print 'rfj rj2 should have done all combinations'

        #third comb rfj jet and H1
            if(recojet_rfj_vector_combto3[2].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3[2].Angle(tempH2P4.Vect()) and recojet_rfj_vector_combto3[2].Angle(tempH1P4.Vect())<recojet_rfj_vector_combto3[2].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[229].Fill(recojet_rfj_vector_combto3[2].P()/tempH1P4.P(),weight)
        #third comb rfj jet and H2
            elif(recojet_rfj_vector_combto3[2].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3[2].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3[2].Angle(tempH2P4.Vect())<recojet_rfj_vector_combto3[2].Angle(tempZP4.Vect())):
                hist_vec_reco_1D_parton[229].Fill(recojet_rfj_vector_combto3[2].P()/tempH2P4.P(),weight)
        #third comb rfj jet and Z
            elif(recojet_rfj_vector_combto3[2].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3[2].Angle(tempH1P4.Vect()) and recojet_rfj_vector_combto3[2].Angle(tempZP4.Vect())<recojet_rfj_vector_combto3[2].Angle(tempH2P4.Vect())):
                hist_vec_reco_1D_parton[229].Fill(recojet_rfj_vector_combto3[2].P()/tempZP4.P(),weight)
            else:
                print 'rfj rj3 should have done all combinations'


 

        #now histos of H1 and combination of rfj vectors
            if(tempH1P4.Angle(recojet_rfj_vector_combto3[0].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3[1].Vect()) and tempH1P4.Angle(recojet_rfj_vector_combto3[0].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3[2].Vect())):
                hist_vec_reco_1D_parton[236].Fill(recojet_rfj_vector_combto3[0].P()/tempH1P4.P(),weight)
                hist_vec_reco_1D_parton[245].Fill(degrees(tempH1P4.Angle(recojet_rfj_vector_combto3[0].Vect())),weight)
            elif(tempH1P4.Angle(recojet_rfj_vector_combto3[1].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3[0].Vect()) and tempH1P4.Angle(recojet_rfj_vector_combto3[1].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3[2].Vect())):
                hist_vec_reco_1D_parton[236].Fill(recojet_rfj_vector_combto3[1].P()/tempH1P4.P(),weight)
                hist_vec_reco_1D_parton[245].Fill(degrees(tempH1P4.Angle(recojet_rfj_vector_combto3[1].Vect())),weight)
            elif(tempH1P4.Angle(recojet_rfj_vector_combto3[2].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3[0].Vect()) and tempH1P4.Angle(recojet_rfj_vector_combto3[2].Vect())<tempH1P4.Angle(recojet_rfj_vector_combto3[1].Vect())):
                hist_vec_reco_1D_parton[236].Fill(recojet_rfj_vector_combto3[2].P()/tempH1P4.P(),weight)
                hist_vec_reco_1D_parton[245].Fill(degrees(tempH1P4.Angle(recojet_rfj_vector_combto3[2].Vect())),weight)
            else:
                print 'sth wrong in H1 to rfj rj comb lineup'

       #now histos of H2 and combination of rfj vectors
            if(tempH2P4.Angle(recojet_rfj_vector_combto3[0].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3[1].Vect()) and tempH2P4.Angle(recojet_rfj_vector_combto3[0].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3[2].Vect())):
                hist_vec_reco_1D_parton[237].Fill(recojet_rfj_vector_combto3[0].P()/tempH2P4.P(),weight)
                hist_vec_reco_1D_parton[246].Fill(degrees(tempH2P4.Angle(recojet_rfj_vector_combto3[0].Vect())),weight)
            elif(tempH2P4.Angle(recojet_rfj_vector_combto3[1].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3[0].Vect()) and tempH2P4.Angle(recojet_rfj_vector_combto3[1].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3[2].Vect())):
                hist_vec_reco_1D_parton[237].Fill(recojet_rfj_vector_combto3[1].P()/tempH2P4.P(),weight)
                hist_vec_reco_1D_parton[246].Fill(degrees(tempH2P4.Angle(recojet_rfj_vector_combto3[1].Vect())),weight)
            elif(tempH2P4.Angle(recojet_rfj_vector_combto3[2].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3[0].Vect()) and tempH2P4.Angle(recojet_rfj_vector_combto3[2].Vect())<tempH2P4.Angle(recojet_rfj_vector_combto3[1].Vect())):
                hist_vec_reco_1D_parton[237].Fill(recojet_rfj_vector_combto3[2].P()/tempH2P4.P(),weight)
                hist_vec_reco_1D_parton[246].Fill(degrees(tempH2P4.Angle(recojet_rfj_vector_combto3[2].Vect())),weight)
            else:
                print 'sth wrong in H2 to rfj rj comb lineup'

      #now histos of Z and combination of rfj vectors
            if(tempZP4.Angle(recojet_rfj_vector_combto3[0].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3[1].Vect()) and tempZP4.Angle(recojet_rfj_vector_combto3[0].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3[2].Vect())):
                hist_vec_reco_1D_parton[238].Fill(recojet_rfj_vector_combto3[0].P()/tempZP4.P(),weight)
                hist_vec_reco_1D_parton[247].Fill(degrees(tempZP4.Angle(recojet_rfj_vector_combto3[0].Vect())),weight)
            elif(tempZP4.Angle(recojet_rfj_vector_combto3[1].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3[0].Vect()) and tempZP4.Angle(recojet_rfj_vector_combto3[1].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3[2].Vect())):
                hist_vec_reco_1D_parton[238].Fill(recojet_rfj_vector_combto3[1].P()/tempZP4.P(),weight)
                hist_vec_reco_1D_parton[247].Fill(degrees(tempZP4.Angle(recojet_rfj_vector_combto3[1].Vect())),weight)
            elif(tempZP4.Angle(recojet_rfj_vector_combto3[2].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3[0].Vect()) and tempZP4.Angle(recojet_rfj_vector_combto3[2].Vect())<tempZP4.Angle(recojet_rfj_vector_combto3[1].Vect())):
                hist_vec_reco_1D_parton[238].Fill(recojet_rfj_vector_combto3[2].P()/tempZP4.P(),weight)
                hist_vec_reco_1D_parton[247].Fill(degrees(tempZP4.Angle(recojet_rfj_vector_combto3[2].Vect())),weight)
            else:
                print 'sth wrong in Z to rfj rj comb lineup'


            if len(genjet_vector)>2: 
                if(tempH1P4.E()>tempH2P4.E()):
                    hist_vec_reco_1D_parton[101].Fill(degrees(genjet_vector[0].Angle(tempH1P4.Vect())),weight)
                    hist_vec_reco_1D_parton[102].Fill(degrees(genjet_vector[0].Angle(tempH2P4.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[107].Fill(degrees(genjet_vector[1].Angle(tempH1P4.Vect())),weight)
                    hist_vec_reco_1D_parton[108].Fill(degrees(genjet_vector[1].Angle(tempH2P4.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[83].Fill(degrees(genjet_vector[2].Angle(tempH1P4.Vect())),weight)
                    hist_vec_reco_1D_parton[84].Fill(degrees(genjet_vector[2].Angle(tempH2P4.Vect())),weight)
                else:
                    hist_vec_reco_1D_parton[101].Fill(degrees(genjet_vector[0].Angle(tempH2P4.Vect())),weight)
                    hist_vec_reco_1D_parton[102].Fill(degrees(genjet_vector[0].Angle(tempH1P4.Vect())),weight)

                    hist_vec_reco_1D_parton[104].Fill(degrees(genjet_vector[1].Angle(tempH2P4.Vect())),weight)
                    hist_vec_reco_1D_parton[105].Fill(degrees(genjet_vector[1].Angle(tempH1P4.Vect())),weight)
                    
                    hist_vec_reco_1D_parton[107].Fill(degrees(genjet_vector[2].Angle(tempH2P4.Vect())),weight)
                    hist_vec_reco_1D_parton[108].Fill(degrees(genjet_vector[2].Angle(tempH1P4.Vect())),weight)
                    
                hist_vec_reco_1D_parton[103].Fill(degrees(genjet_vector[1].Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[106].Fill(degrees(genjet_vector[1].Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[109].Fill(degrees(genjet_vector[2].Angle(tempZP4.Vect())),weight)
                
                hist_vec_reco_1D_parton[110].Fill(genjet_vector[0].M(),weight)
                hist_vec_reco_1D_parton[111].Fill(genjet_vector[1].M(),weight)
                hist_vec_reco_1D_parton[112].Fill(genjet_vector[2].M(),weight)
            
                hist_vec_reco_1D_parton[129].Fill(genjet_vector[0].E(),weight)
                hist_vec_reco_1D_parton[130].Fill(genjet_vector[1].E(),weight)
                hist_vec_reco_1D_parton[131].Fill(genjet_vector[2].E(),weight)
            
                hist_vec_reco_1D_parton[151].Fill(genjet_vector_combto3[0].M(),weight)
                hist_vec_reco_1D_parton[152].Fill(genjet_vector_combto3[1].M(),weight)
                hist_vec_reco_1D_parton[153].Fill(genjet_vector_combto3[2].M(),weight)

            if len(genjet_vector)>3:
                hist_vec_reco_1D_parton[123].Fill(genjet_vector[3].M(),weight)
                hist_vec_reco_1D_parton[125].Fill(degrees(min(genjet_vector[3].Angle(tempH1P4.Vect()),genjet_vector[3].Angle(tempH2P4.Vect()))),weight)
                hist_vec_reco_1D_parton[126].Fill(degrees(genjet_vector[3].Angle(tempZP4.Vect())),weight)
                hist_vec_reco_1D_parton[132].Fill(genjet_vector[3].E(),weight)
                if len(genjet_vector)>4:
                    hist_vec_reco_1D_parton[137].Fill(degrees(min(genjet_vector[4].Angle(tempH1P4.Vect()),genjet_vector[4].Angle(tempH2P4.Vect()))),weight)
                    hist_vec_reco_1D_parton[138].Fill(degrees(genjet_vector[4].Angle(tempZP4.Vect())),weight)
                    if len(genjet_vector)>5:
                        hist_vec_reco_1D_parton[141].Fill(degrees(min(genjet_vector[5].Angle(tempH1P4.Vect()),genjet_vector[5].Angle(tempH2P4.Vect()))),weight)
                        hist_vec_reco_1D_parton[142].Fill(degrees(genjet_vector[5].Angle(tempZP4.Vect())),weight)

    if len(hist_vec_reco_1D_parton)>0:
        print 'total events after all running signal histos', hist_vec_reco_1D_parton[0].Integral(0,hist_vec_reco_1D_parton[0].GetNbinsX()+1)," ",hist_vec_reco_1D_parton[110].Integral(0,hist_vec_reco_1D_parton[110].GetNbinsX()+1),num_entry,hist_vec_reco_1D_parton[0].GetEntries(),num_count, num_total_exception

    return None


    

def process_event(i_final_histo_name_,i_input_file_name_,i_xsec_,i_lumi_,i_fillPartonHistos,i_ishhzfile):
    print "at start of process event"         
    input_file_=root.TFile.Open(i_input_file_name_)
    #input_file2_=root.TFile.Open(i_input_file_name2_)
    lumi=i_lumi_
    xsec_=i_xsec_
    print 'process event',i_xsec_,i_lumi_,i_fillPartonHistos,i_ishhzfile
          
                            
    file_histogram = root.TFile(i_final_histo_name_, "RECREATE")

    _mytree = TTree('MVATrainingVariables', 'MVATrainingVariables')


    hist_vec_HHZ_parton_list=[]
 
    if i_fillPartonHistos:
        n_bins_high=200;
        
        lim_energy_low=0.;
        lim_energy_high=3050.;
        
        h_sqrtS_e1_e2_effective = TH1F("h_sqrtS_e1_e2_effective","", n_bins_high, lim_energy_low,lim_energy_high);
        h_sqrtS_e1_e2_effective_HHZ_bbbbqq = TH1F("h_sqrtS_e1_e2_effective_HHZ_bbbbqq","", n_bins_high, lim_energy_low,lim_energy_high);
        h_sqrtS_HHZ = TH1F("h_sqrtS_HHZ","", n_bins_high, lim_energy_low,lim_energy_high);
        h_H1_H2_E_sum = TH1F("h_H1_H2_E_sum","", n_bins_high, lim_energy_low,lim_energy_high);
        h_H1_E = TH1F("h_H1_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_H2_E = TH1F("h_H2_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_Z_E = TH1F("h_Z_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        
        h_H1_P = TH1F("h_H1_P","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_H2_P = TH1F("h_H2_P","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_Z_P = TH1F("h_Z_P","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        
        h_sqrtS_gj_sum_all = TH1F("h_sqrtS_gj_sum_all","", n_bins_high, lim_energy_low,lim_energy_high);
        h_sqrtS_gj_sum_HHZ_bbbbqq = TH1F("h_sqrtS_gj_sum_HHZ_bbbbqq ","", n_bins_high, lim_energy_low,lim_energy_high);
        h_sqrtS_gj_sum_plus_InvVec_all = TH1F("h_sqrtS_gj_sum_plus_InvVec_all","", n_bins_high, lim_energy_low,lim_energy_high);
        h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq = TH1F("h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq ","", n_bins_high, lim_energy_low,lim_energy_high);
        
        lim_mass_low=0;
        lim_mass_high=500;
        
        h_H1_H2_mass = TH1F("h_H1_H2_mass","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        
        lim_dalpha_low=0.;
        lim_dalpha_high=180.;
        
        h_theta_H1 = TH1F("h_theta_H1","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_theta_H2 = TH1F("h_theta_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_theta_Z = TH1F("h_theta_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        
        h_dalpha_H1_H2_comb_vs_Z = TH1F("h_dalpha_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dphi_H1_H2_comb_vs_Z = TH1F("h_dphi_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dtheta_H1_H2_comb_vs_Z = TH1F("h_dtheta_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        
        h_dalpha_H1_H2 = TH1F("h_dalpha_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dphi_H1_H2 = TH1F("h_dphi_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dtheta_H1_H2 = TH1F("h_dtheta_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        
        h_dalpha_H_Z_min = TH1F("h_dalpha_H_Z_min","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dphi_H_Z_min = TH1F("h_dphi_H_Z_min","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dtheta_H_Z_min = TH1F("h_dtheta_H_Z_min","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dalpha_H_Z_max = TH1F("h_dalpha_H_Z_max","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dphi_H_Z_max = TH1F("h_dphi_H_Z_max","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dtheta_H_Z_max = TH1F("h_dtheta_H_Z_max","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_Esum_Z_H_minAngle = TH1F("Esum_Z_H_minAngle","", n_bins_high, lim_energy_low,lim_energy_high);
        h_delta_Esum_Z_H_minAngle_min_E_H_Z_maxAngle = TH1F("h_delta_Esum_Z_H_minAngle_min_E_H_Z_maxAngle","", n_bins_high, -0.5*lim_energy_high,0.5*lim_energy_high);
        
        lim_dalpha_qqbar_low=0.;
        lim_dalpha_qqbar_high=45.;
        
        h_dalpha_H1_bbbar = TH1F("h_dalpha_H1_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dphi_H1_bbbar = TH1F("h_dphi_H1_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dtheta_H1_bbbar = TH1F("h_dtheta_H1_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        
        h_dalpha_H2_bbbar = TH1F("h_dalpha_H2_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dphi_H2_bbbar = TH1F("h_dphi_H2_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dtheta_H2_bbbar = TH1F("h_dtheta_H2_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        
        h_dalpha_H1_qqbar = TH1F("h_dalpha_H1_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dphi_H1_qqbar = TH1F("h_dphi_H1_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dtheta_H1_qqbar = TH1F("h_dtheta_H1_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        
        h_dalpha_H2_qqbar = TH1F("h_dalpha_H2_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dphi_H2_qqbar = TH1F("h_dphi_H2_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dtheta_H2_qqbar = TH1F("h_dtheta_H2_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        
        h_dalpha_Z_qqbar = TH1F("h_dalpha_Z_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dphi_Z_qqbar = TH1F("h_dphi_Z_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dtheta_Z_qqbar = TH1F("h_dtheta_Z_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        
        
        lim_dalpha_qqbar_allH_low=0.;
        lim_dalpha_qqbar_allH_high=100.;
        
        h_dalpha_max_allH_qqbar = TH1F("h_dalpha_max_allH_qqbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
        h_dphi_max_allH_qqbar = TH1F("h_dphi_max_allH_qqbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
        h_dtheta_max_allH_qqbar = TH1F("h_dtheta_max_allH_qqbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
        
        h_dalpha_max_allH_bbbar = TH1F("h_dalpha_max_allH_bbbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
        h_dphi_max_allH_bbbar = TH1F("h_dphi_max_allH_bbbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
        h_dtheta_max_allH_bbbar = TH1F("h_dtheta_max_allH_bbbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
        
        h_dalpha_allH_diboson_qqbar = TH1F("h_dalpha_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dphi_allH_diboson_qqbar = TH1F("h_dphi_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        h_dtheta_allH_diboson_qqbar = TH1F("h_dtheta_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
        
        lim_H1_E_over_allH_E_low=0.5;
        lim_H1_E_over_allH_E_high=1.0;
        
        h_E_H1_over_E_allH = TH1F("h_E_H1_over_E_allH","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
        h_E_H1_over_E_allH_bbbar = TH1F("h_E_H1_over_E_allH_bbbar","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
        h_E_H1_over_E_allH_qqbar = TH1F("h_E_H1_over_E_allH_qqbar","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
        
        lim_d_ij_qqbar_low=0.;
        lim_d_ij_qqbar_high=5.0;
        
        h_d_ij_H1_bbbar = TH1F("h_d_ij_H1_bbbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        h_d_ij_H2_bbbar = TH1F("h_d_ij_H2_bbbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        
        h_d_ij_H1_H2 = TH1F("h_d_ij_H1_H2","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        h_d_ij_H1_Z = TH1F("h_d_ij_H1_Z","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        h_d_ij_H2_Z = TH1F("h_d_ij_H2_Z","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        h_d_ij_H1_qqbar = TH1F("h_d_ij_H1_qqbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        h_d_ij_H2_qqbar = TH1F("h_d_ij_H2_qqbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        h_d_ij_Z_qqbar = TH1F("h_d_ij_Z_qqbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        
        h_dalpha_min_V1_V2 = TH1F("h_dalpha_min_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dphi_min_V1_V2 = TH1F("h_dphi_min_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dtheta_min_V1_V2 = TH1F("h_dtheta_min_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        
        h_dalpha_max_V1_V2 = TH1F("h_dalpha_max_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dphi_max_V1_V2 = TH1F("h_dphi_max_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
        h_dtheta_max_V1_V2 = TH1F("h_dtheta_max_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

        h_d_ij_min_V1_V2 = TH1F("h_d_ij_min_V1_V2","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        h_d_ij_max_V1_V2 = TH1F("h_d_ij_max_V1_V2","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
        
        lim_angle_low=0
        lim_angle_high=180
        
        h_dalpha_H1_rj_HHZ_bbbbqq = TH1F("h_dalpha_H1_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dphi_H1_rj_HHZ_bbbbqq = TH1F("h_dphi_H1_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dtheta_H1_rj_HHZ_bbbbqq = TH1F("h_dtheta_H1_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_H2_rj_HHZ_bbbbqq = TH1F("h_dalpha_H2_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dphi_H2_rj_HHZ_bbbbqq = TH1F("h_dphi_H2_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dtheta_H2_rj_HHZ_bbbbqq = TH1F("h_dtheta_H2_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_Z_rj_HHZ_bbbbqq = TH1F("h_dalpha_Z_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dphi_Z_rj_HHZ_bbbbqq = TH1F("h_dphi_Z_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dtheta_Z_rj_HHZ_bbbbqq = TH1F("h_dtheta_Z_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_H1_gj_HHZ_bbbbqq = TH1F("h_dalpha_H1_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dphi_H1_gj_HHZ_bbbbqq = TH1F("h_dphi_H1_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dtheta_H1_gj_HHZ_bbbbqq = TH1F("h_dtheta_H1_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_H2_gj_HHZ_bbbbqq = TH1F("h_dalpha_H2_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dphi_H2_gj_HHZ_bbbbqq = TH1F("h_dphi_H2_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dtheta_H2_gj_HHZ_bbbbqq = TH1F("h_dtheta_H2_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_Z_gj_HHZ_bbbbqq = TH1F("h_dalpha_Z_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dphi_Z_gj_HHZ_bbbbqq = TH1F("h_dphi_Z_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dtheta_Z_gj_HHZ_bbbbqq = TH1F("h_dtheta_Z_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        lim_mass_low=0
        lim_mass_high=250
        h_mass_Z_rj_HHZ_bbbbqq = TH1F("h_mass_Z_rj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_H1_rj_HHZ_bbbbqq = TH1F("h_mass_H1_rj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_H2_rj_HHZ_bbbbqq = TH1F("h_mass_H2_rj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        
        h_dalpha_rj1_H1_HHZ_bbbbqq = TH1F("h_dalpha_rj1_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj1_H2_HHZ_bbbbqq = TH1F("h_dalpha_rj1_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj1_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj1_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_rj2_H1_HHZ_bbbbqq = TH1F("h_dalpha_rj2_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj2_H2_HHZ_bbbbqq = TH1F("h_dalpha_rj2_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj2_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj2_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_rj3_H1_HHZ_bbbbqq = TH1F("h_dalpha_rj3_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj3_H2_HHZ_bbbbqq = TH1F("h_dalpha_rj3_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj3_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj3_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_mass_rj1_HHZ_bbbbqq = TH1F("h_mass_rj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_rj2_HHZ_bbbbqq = TH1F("h_mass_rj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_rj3_HHZ_bbbbqq = TH1F("h_mass_rj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        
        h_mass_Z_gj_HHZ_bbbbqq = TH1F("h_mass_Z_gj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_H1_gj_HHZ_bbbbqq = TH1F("h_mass_H1_gj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_H2_gj_HHZ_bbbbqq = TH1F("h_mass_H2_gj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        
        h_dalpha_gj1_H1_HHZ_bbbbqq = TH1F("h_dalpha_gj1_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_gj1_H2_HHZ_bbbbqq = TH1F("h_dalpha_gj1_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_gj1_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj1_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_gj2_H1_HHZ_bbbbqq = TH1F("h_dalpha_gj2_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_gj2_H2_HHZ_bbbbqq = TH1F("h_dalpha_gj2_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_gj2_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj2_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_gj3_H1_HHZ_bbbbqq = TH1F("h_dalpha_gj3_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_gj3_H2_HHZ_bbbbqq = TH1F("h_dalpha_gj3_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_gj3_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj3_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_mass_gj1_HHZ_bbbbqq = TH1F("h_mass_gj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_gj2_HHZ_bbbbqq = TH1F("h_mass_gj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_gj3_HHZ_bbbbqq = TH1F("h_mass_gj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        
        h_mass_gj4_HHZ_bbbbqq = TH1F("h_mass_gj4_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_rj4_HHZ_bbbbqq = TH1F("h_mass_rj4_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);

        h_dalpha_gj4_Hs_HHZ_bbbbqq = TH1F("h_dalpha_gj4_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_gj4_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj4_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

        h_dalpha_rj4_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rj4_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj4_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj4_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_gj1_E_bbbbqq = TH1F("h_gj1_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_gj2_E_bbbbqq = TH1F("h_gj2_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_gj3_E_bbbbqq = TH1F("h_gj3_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_gj4_E_bbbbqq = TH1F("h_gj4_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        
        h_rj1_E_bbbbqq = TH1F("h_rj1_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_rj2_E_bbbbqq = TH1F("h_rj2_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_rj3_E_bbbbqq = TH1F("h_rj3_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        h_rj4_E_bbbbqq = TH1F("h_rj4_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
        
        h_dalpha_gj5_Hs_HHZ_bbbbqq = TH1F("h_dalpha_gj5_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_gj5_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj5_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_rj5_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rj5_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj5_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj5_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_gj6_Hs_HHZ_bbbbqq = TH1F("h_dalpha_gj6_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_gj6_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj6_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_rj6_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rj6_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj6_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj6_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_rj1_gjs_bbbbqq = TH1F("h_dalpha_rj1_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj2_gjs_bbbbqq = TH1F("h_dalpha_rj2_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj3_gjs_bbbbqq = TH1F("h_dalpha_rj3_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj4_gjs_bbbbqq = TH1F("h_dalpha_rj4_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj5_gjs_bbbbqq = TH1F("h_dalpha_rj5_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rj6_gjs_bbbbqq = TH1F("h_dalpha_rj6_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_mass_comb_gj1_HHZ_bbbbqq = TH1F("h_mass_comb_gj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_comb_gj2_HHZ_bbbbqq = TH1F("h_mass_comb_gj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_comb_gj3_HHZ_bbbbqq = TH1F("h_mass_comb_gj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        
        h_mass_comb_rj1_HHZ_bbbbqq = TH1F("h_mass_comb_rj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_comb_rj2_HHZ_bbbbqq = TH1F("h_mass_comb_rj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_comb_rj3_HHZ_bbbbqq = TH1F("h_mass_comb_rj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        
        h_dalpha_comb_rj1_gjs_bbbbqq = TH1F("h_dalpha_comb_rj1_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rj2_gjs_bbbbqq = TH1F("h_dalpha_comb_rj2_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rj3_gjs_bbbbqq = TH1F("h_dalpha_comb_rj3_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_comb_rj1_Hs_bbbbqq = TH1F("h_dalpha_comb_rj1_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rj1_Z_bbbbqq = TH1F("h_dalpha_comb_rj1_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rj2_Hs_bbbbqq = TH1F("h_dalpha_comb_rj2_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rj2_Z_bbbbqq = TH1F("h_dalpha_comb_rj2_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rj3_Hs_bbbbqq = TH1F("h_dalpha_comb_rj3_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rj3_Z_bbbbqq = TH1F("h_dalpha_comb_rj3_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_mass_comb_rfj_rj1_HHZ_bbbbqq = TH1F("h_mass_comb_rfj_rj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_comb_rfj_rj2_HHZ_bbbbqq = TH1F("h_mass_comb_rfj_rj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_comb_rfj_rj3_HHZ_bbbbqq = TH1F("h_mass_comb_rfj_rj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        
        h_dalpha_comb_rfj_rj1_gjs_bbbbqq = TH1F("h_dalpha_comb_rfj_rj1_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rfj_rj2_gjs_bbbbqq = TH1F("h_dalpha_comb_rfj_rj2_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rfj_rj3_gjs_bbbbqq = TH1F("h_dalpha_comb_rfj_rj3_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_comb_rfj_rj1_Hs_bbbbqq = TH1F("h_dalpha_comb_rfj_rj1_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rfj_rj1_Z_bbbbqq = TH1F("h_dalpha_comb_rfj_rj1_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rfj_rj2_Hs_bbbbqq = TH1F("h_dalpha_comb_rfj_rj2_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rfj_rj2_Z_bbbbqq = TH1F("h_dalpha_comb_rfj_rj2_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rfj_rj3_Hs_bbbbqq = TH1F("h_dalpha_comb_rfj_rj3_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_rfj_rj3_Z_bbbbqq = TH1F("h_dalpha_comb_rfj_rj3_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        lim_BTagSum_low=0
        lim_BTagSum_high=2
        
        h_BTagSum_comb_rfj_rj1_HHZ_bbbbqq = TH1F("h_BTagSum_comb_rfj_rj1_HHZ_bbbbqq","", n_bins_high, lim_BTagSum_low,lim_BTagSum_high);
        h_BTagSum_comb_rfj_rj2_HHZ_bbbbqq = TH1F("h_BTagSum_comb_rfj_rj2_HHZ_bbbbqq","", n_bins_high, lim_BTagSum_low,lim_BTagSum_high);
        h_BTagSum_comb_rfj_rj3_HHZ_bbbbqq = TH1F("h_BTagSum_comb_rfj_rj3_HHZ_bbbbqq","", n_bins_high, lim_BTagSum_low,lim_BTagSum_high);
        
        lim_BTagMax_low=0
        lim_BTagMax_high=1
        
        h_BTagMax_comb_rfj_rj1_HHZ_bbbbqq = TH1F("h_BTagMax_comb_rfj_rj1_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        h_BTagMax_comb_rfj_rj2_HHZ_bbbbqq = TH1F("h_BTagMax_comb_rfj_rj2_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        h_BTagMax_comb_rfj_rj3_HHZ_bbbbqq = TH1F("h_BTagMax_comb_rfj_rj3_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        
        
        h_BTagMax_rfj1_BTag_HHZ_bbbbqq = TH1F("h_BTagMax_rfj1_BTag_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        h_BTagMax_rfj2_BTag_HHZ_bbbbqq = TH1F("h_BTagMax_rfj2_BTag_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        h_BTagMax_rfj3_BTag_HHZ_bbbbqq = TH1F("h_BTagMax_rfj3_BTag_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        h_BTagMax_rfj4_BTag_HHZ_bbbbqq = TH1F("h_BTagMax_rfj4_BTag_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        h_BTagMax_rfj5_BTag_HHZ_bbbbqq = TH1F("h_BTagMax_rfj5_BTag_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        h_BTagMax_rfj6_BTag_HHZ_bbbbqq = TH1F("h_BTagMax_rfj6_BTag_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        
        h_dalpha_rfj_low_BTag_highMass_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rfj_low_BTag_highMass_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rfj_low_BTag_highMass_Z_HHZ_bbbbqq = TH1F("h_dalpha_rfj_low_BTag_highMass_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rfj_low_BTag_highMass_b_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rfj_low_BTag_highMass_b_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rfj_low_BTag_highMass_q_Z_HHZ_bbbbqq = TH1F("h_dalpha_rfj_low_BTag_highMass_q_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_rfj_high_BTag_highMass_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rfj_high_BTag_highMass_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rfj_high_BTag_highMass_Z_HHZ_bbbbqq = TH1F("h_dalpha_rfj_high_BTag_highMass_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rfj_high_BTag_highMass_b_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rfj_high_BTag_highMass_b_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rfj_high_BTag_highMass_q_Z_HHZ_bbbbqq = TH1F("h_dalpha_rfj_high_BTag_highMass_q_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_rfj_low_BTag_highMass_gj_HHZ_bbbbqq = TH1F("h_dalpha_rfj_low_BTag_highMass_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_rfj_high_BTag_highMass_gj_HHZ_bbbbqq = TH1F("h_dalpha_rfj_high_BTag_highMass_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        
        n_bins_high_extreme=3000
        lim_mass_low_extreme=0
        lim_mass_high_extreme=1500
        
        h_mass_comb_rj1_HHZ_bbbbqq_forRMS90 = TH1F("h_mass_comb_rj1_HHZ_bbbbqq_forRMS90","", n_bins_high_extreme, lim_mass_low_extreme,lim_mass_high_extreme);
        h_mass_comb_rj2_HHZ_bbbbqq_forRMS90 = TH1F("h_mass_comb_rj2_HHZ_bbbbqq_forRMS90","", n_bins_high_extreme, lim_mass_low_extreme,lim_mass_high_extreme);
        h_mass_comb_rj3_HHZ_bbbbqq_forRMS90 = TH1F("h_mass_comb_rj3_HHZ_bbbbqq_forRMS90","", n_bins_high_extreme, lim_mass_low_extreme,lim_mass_high_extreme);
        
        h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90 = TH1F("h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90","", n_bins_high_extreme, lim_mass_low_extreme,lim_mass_high_extreme);
        h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90 = TH1F("h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90","", n_bins_high_extreme, lim_mass_low_extreme,lim_mass_high_extreme);
        h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90 = TH1F("h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90","", n_bins_high_extreme, lim_mass_low_extreme,lim_mass_high_extreme);
        
        h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90 = TH1F("h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90","", n_bins_high_extreme, lim_mass_low_extreme,lim_mass_high_extreme);
        h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90 = TH1F("h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90","", n_bins_high_extreme, lim_mass_low_extreme,lim_mass_high_extreme);
        h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90 = TH1F("h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90","", n_bins_high_extreme, lim_mass_low_extreme,lim_mass_high_extreme);
        
        h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq = TH1F("h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq = TH1F("h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq = TH1F("h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
        
        h_dalpha_comb_wBTag_rfj_rj1_gjs_bbbbqq = TH1F("h_dalpha_comb_wBTag_rfj_rj1_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_wBTag_rfj_rj2_gjs_bbbbqq = TH1F("h_dalpha_comb_wBTag_rfj_rj2_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_wBTag_rfj_rj3_gjs_bbbbqq = TH1F("h_dalpha_comb_wBTag_rfj_rj3_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_comb_wBTag_rfj_rj1_Hs_bbbbqq = TH1F("h_dalpha_comb_wBTag_rfj_rj1_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_wBTag_rfj_rj1_Z_bbbbqq = TH1F("h_dalpha_comb_wBTag_rfj_rj1_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_wBTag_rfj_rj2_Hs_bbbbqq = TH1F("h_dalpha_comb_wBTag_rfj_rj2_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_wBTag_rfj_rj2_Z_bbbbqq = TH1F("h_dalpha_comb_wBTag_rfj_rj2_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_wBTag_rfj_rj3_Hs_bbbbqq = TH1F("h_dalpha_comb_wBTag_rfj_rj3_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_comb_wBTag_rfj_rj3_Z_bbbbqq = TH1F("h_dalpha_comb_wBTag_rfj_rj3_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_BTagMax_comb_wBTag_rfj_rj1_HHZ_bbbbqq = TH1F("h_BTagMax_comb_wBTag_rfj_rj1_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        h_BTagMax_comb_wBTag_rfj_rj2_HHZ_bbbbqq = TH1F("h_BTagMax_comb_wBTag_rfj_rj2_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        h_BTagMax_comb_wBTag_rfj_rj3_HHZ_bbbbqq = TH1F("h_BTagMax_comb_wBTag_rfj_rj3_HHZ_bbbbqq","", n_bins_high, lim_BTagMax_low,lim_BTagMax_high);
        
        lim_P_part_over_reco_low=0
        lim_P_part_over_reco_high=4.0
        h_comb_rj1_matchPart_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_comb_rj1_matchPart_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_comb_rj2_matchPart_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_comb_rj2_matchPart_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_comb_rj3_matchPart_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_comb_rj3_matchPart_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        
        h_comb_rfj_rj1_matchPart_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_comb_rfj_rj1_matchPart_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_comb_rfj_rj2_matchPart_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_comb_rfj_rj2_matchPart_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_comb_rfj_rj3_matchPart_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_comb_rfj_rj3_matchPart_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        
        h_comb_wBTag_rfj_rj1_matchPart_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_comb_wBTag_rfj_rj1_matchPart_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_comb_wBTag_rfj_rj2_matchPart_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_comb_wBTag_rfj_rj2_matchPart_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_comb_wBTag_rfj_rj3_matchPart_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_comb_wBTag_rfj_rj3_matchPart_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        
        h_H1_comb_rj_match_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_H1_comb_rj_match_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_H2_comb_rj_match_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_H2_comb_rj_match_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_Z_comb_rj_match_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_Z_comb_rj_match_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        
        h_H1_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_H1_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_H2_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_H2_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_Z_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_Z_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        
        h_H1_comb_wBTag_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_H1_comb_wBTag_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_H2_comb_wBTag_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_H2_comb_wBTag_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        h_Z_comb_wBTag_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq = TH1F("h_Z_comb_wBTag_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq","", n_bins_high, lim_P_part_over_reco_low,lim_P_part_over_reco_high);
        
        h_dalpha_H1_comb_rj_match_HHZ_bbbbqq = TH1F("h_dalpha_H1_comb_rj_match_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_H2_comb_rj_match_HHZ_bbbbqq = TH1F("h_dalpha_H2_comb_rj_match_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_Z_comb_rj_match_HHZ_bbbbqq = TH1F("h_dalpha_Z_comb_rj_match_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_H1_comb_rfj_rj_match_HHZ_bbbbqq = TH1F("h_dalpha_H1_comb_rfj_rj_match_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_H2_comb_rfj_rj_match_HHZ_bbbbqq = TH1F("h_dalpha_H2_comb_rfj_rj_match_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq = TH1F("h_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        
        h_dalpha_H1_comb_wBTag_rfj_rj_match_HHZ_bbbbqq = TH1F("h_dalpha_H1_comb_wBTag_rfj_rj_match_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_H2_comb_wBTag_rfj_rj_match_HHZ_bbbbqq = TH1F("h_dalpha_H2_comb_wBTag_rfj_rj_match_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        h_dalpha_Z_comb_wBTag_rfj_rj_match_HHZ_bbbbqq = TH1F("h_dalpha_Z_comb_wBTag_rfj_rj_match_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
        

        
        hist_vec_HHZ_parton_list.append(h_sqrtS_e1_e2_effective)
        hist_vec_HHZ_parton_list.append(h_H1_H2_E_sum)
        hist_vec_HHZ_parton_list.append(h_H1_E)
        hist_vec_HHZ_parton_list.append(h_H2_E)
        hist_vec_HHZ_parton_list.append(h_Z_E)
    #4 done here
        hist_vec_HHZ_parton_list.append(h_H1_H2_mass)
        hist_vec_HHZ_parton_list.append(h_dalpha_H1_H2_comb_vs_Z)
        hist_vec_HHZ_parton_list.append(h_dphi_H1_H2_comb_vs_Z)
        hist_vec_HHZ_parton_list.append(h_dtheta_H1_H2_comb_vs_Z)
        hist_vec_HHZ_parton_list.append(h_dalpha_H1_H2)
    #9 done here
        hist_vec_HHZ_parton_list.append(h_dphi_H1_H2)
        hist_vec_HHZ_parton_list.append(h_dtheta_H1_H2)
        hist_vec_HHZ_parton_list.append(h_dalpha_H1_bbbar)
        hist_vec_HHZ_parton_list.append(h_dphi_H1_bbbar)
        hist_vec_HHZ_parton_list.append(h_dtheta_H1_bbbar)
    #14 done here
        hist_vec_HHZ_parton_list.append(h_dalpha_H2_bbbar)
        hist_vec_HHZ_parton_list.append(h_dphi_H2_bbbar)
        hist_vec_HHZ_parton_list.append(h_dtheta_H2_bbbar)
        hist_vec_HHZ_parton_list.append(h_dalpha_H1_qqbar)
        hist_vec_HHZ_parton_list.append(h_dphi_H1_qqbar)
    #19 done here
        hist_vec_HHZ_parton_list.append(h_dtheta_H1_qqbar)
        hist_vec_HHZ_parton_list.append(h_dalpha_H2_qqbar)
        hist_vec_HHZ_parton_list.append(h_dphi_H2_qqbar)
        hist_vec_HHZ_parton_list.append(h_dtheta_H2_qqbar)
        hist_vec_HHZ_parton_list.append(h_dalpha_Z_qqbar)
    #24 done here
        hist_vec_HHZ_parton_list.append(h_dphi_Z_qqbar)
        hist_vec_HHZ_parton_list.append(h_dtheta_Z_qqbar)
        hist_vec_HHZ_parton_list.append(h_dalpha_max_allH_qqbar)
        hist_vec_HHZ_parton_list.append(h_dphi_max_allH_qqbar)
        hist_vec_HHZ_parton_list.append(h_dtheta_max_allH_qqbar)
    #29 done here
        hist_vec_HHZ_parton_list.append(h_dalpha_max_allH_bbbar)
        hist_vec_HHZ_parton_list.append(h_dphi_max_allH_bbbar)
        hist_vec_HHZ_parton_list.append(h_dtheta_max_allH_bbbar)
        hist_vec_HHZ_parton_list.append(h_dalpha_allH_diboson_qqbar)
        hist_vec_HHZ_parton_list.append(h_dphi_allH_diboson_qqbar)
    #34 done here
        hist_vec_HHZ_parton_list.append(h_dtheta_allH_diboson_qqbar)
        hist_vec_HHZ_parton_list.append(h_E_H1_over_E_allH)
        hist_vec_HHZ_parton_list.append(h_E_H1_over_E_allH_bbbar)
        hist_vec_HHZ_parton_list.append(h_E_H1_over_E_allH_qqbar)
        hist_vec_HHZ_parton_list.append(h_d_ij_H1_H2)
    #39 done up to here
        hist_vec_HHZ_parton_list.append(h_d_ij_H1_qqbar)
        hist_vec_HHZ_parton_list.append(h_d_ij_H2_qqbar)
        hist_vec_HHZ_parton_list.append(h_sqrtS_HHZ)
        hist_vec_HHZ_parton_list.append(h_dalpha_H_Z_min)
        hist_vec_HHZ_parton_list.append(h_dphi_H_Z_min)
    #44 done up to here
        hist_vec_HHZ_parton_list.append(h_dtheta_H_Z_min)
        hist_vec_HHZ_parton_list.append(h_dalpha_H_Z_max)
        hist_vec_HHZ_parton_list.append(h_dphi_H_Z_max)
        hist_vec_HHZ_parton_list.append(h_dtheta_H_Z_max)
        hist_vec_HHZ_parton_list.append(h_Esum_Z_H_minAngle)
    #49 done up to here --> dphi and dtheta for maximally angularly separated, but still check if it is more on phi or theta
        hist_vec_HHZ_parton_list.append(h_delta_Esum_Z_H_minAngle_min_E_H_Z_maxAngle)
        hist_vec_HHZ_parton_list.append(h_dalpha_min_V1_V2)
        hist_vec_HHZ_parton_list.append(h_dphi_min_V1_V2)
        hist_vec_HHZ_parton_list.append(h_dtheta_min_V1_V2)
        hist_vec_HHZ_parton_list.append(h_dalpha_max_V1_V2)
    #54 done up to now
        hist_vec_HHZ_parton_list.append(h_dphi_max_V1_V2)
        hist_vec_HHZ_parton_list.append(h_dtheta_max_V1_V2)
        hist_vec_HHZ_parton_list.append(h_d_ij_min_V1_V2)
        hist_vec_HHZ_parton_list.append(h_d_ij_max_V1_V2)
        hist_vec_HHZ_parton_list.append(h_d_ij_Z_qqbar)
    #59 done up to now
        hist_vec_HHZ_parton_list.append(h_d_ij_H1_bbbar)
        hist_vec_HHZ_parton_list.append(h_d_ij_H2_bbbar)
        hist_vec_HHZ_parton_list.append(h_theta_H1)
        hist_vec_HHZ_parton_list.append(h_theta_H2)
        hist_vec_HHZ_parton_list.append(h_theta_Z)
    #64 done up to now

        hist_vec_HHZ_parton_list.append(h_dalpha_H1_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dphi_H1_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dtheta_H1_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_H2_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dphi_H2_rj_HHZ_bbbbqq)
    #69 done up to now
        hist_vec_HHZ_parton_list.append(h_dtheta_H2_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_Z_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dphi_Z_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dtheta_Z_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_Z_rj_HHZ_bbbbqq)
    #74 done up to now
        hist_vec_HHZ_parton_list.append(h_mass_H1_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_H2_rj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj1_H1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj1_H2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj1_Z_HHZ_bbbbqq) 
    #79 done up to now
        hist_vec_HHZ_parton_list.append(h_dalpha_rj2_H1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj2_H2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj2_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj3_H1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj3_H2_HHZ_bbbbqq)
    #84 done up to now
        hist_vec_HHZ_parton_list.append(h_dalpha_rj3_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_rj1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_rj2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_rj3_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_H1_gj_HHZ_bbbbqq)
    #89 up to here
        hist_vec_HHZ_parton_list.append(h_dphi_H1_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dtheta_H1_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_H2_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dphi_H2_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dtheta_H2_gj_HHZ_bbbbqq)
    #94 done up to now    
        hist_vec_HHZ_parton_list.append(h_dalpha_Z_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dphi_Z_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dtheta_Z_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_Z_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_H1_gj_HHZ_bbbbqq)
    #99 done up to now   
        hist_vec_HHZ_parton_list.append(h_mass_H2_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj1_H1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj1_H2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj1_Z_HHZ_bbbbqq) 
        hist_vec_HHZ_parton_list.append(h_dalpha_gj2_H1_HHZ_bbbbqq)
    #104 done up to now
        hist_vec_HHZ_parton_list.append(h_dalpha_gj2_H2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj2_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj3_H1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj3_H2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj3_Z_HHZ_bbbbqq)
    #109 done up to now
        hist_vec_HHZ_parton_list.append(h_mass_gj1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_gj2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_gj3_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_sqrtS_gj_sum_all)
        hist_vec_HHZ_parton_list.append(h_sqrtS_gj_sum_HHZ_bbbbqq)
    #114 done until now
        hist_vec_HHZ_parton_list.append(h_sqrtS_gj_sum_plus_InvVec_all)
        hist_vec_HHZ_parton_list.append(h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_sqrtS_e1_e2_effective_HHZ_bbbbqq)

        hist_vec_HHZ_parton_list.append(h_H1_P)
        hist_vec_HHZ_parton_list.append(h_H2_P)
    #119 done until now
        hist_vec_HHZ_parton_list.append(h_Z_P)
        hist_vec_HHZ_parton_list.append(h_d_ij_H1_Z)
        hist_vec_HHZ_parton_list.append(h_d_ij_H2_Z)
        
        hist_vec_HHZ_parton_list.append(h_mass_gj4_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_rj4_HHZ_bbbbqq)
    #124 done until now
        hist_vec_HHZ_parton_list.append(h_dalpha_gj4_Hs_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj4_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj4_Hs_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj4_Z_HHZ_bbbbqq)
        
        hist_vec_HHZ_parton_list.append(h_gj1_E_bbbbqq)
    #129 done until now
        hist_vec_HHZ_parton_list.append(h_gj2_E_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_gj3_E_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_gj4_E_bbbbqq)
        
        hist_vec_HHZ_parton_list.append(h_rj1_E_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_rj2_E_bbbbqq)
    #134 done until now
        hist_vec_HHZ_parton_list.append(h_rj3_E_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_rj4_E_bbbbqq)
        
        hist_vec_HHZ_parton_list.append(h_dalpha_gj5_Hs_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj5_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj5_Hs_HHZ_bbbbqq)
    #139 done
        hist_vec_HHZ_parton_list.append(h_dalpha_rj5_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj6_Hs_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_gj6_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj6_Hs_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj6_Z_HHZ_bbbbqq)
    #144 done 
        hist_vec_HHZ_parton_list.append(h_dalpha_rj1_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj2_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj3_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj4_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rj5_gjs_bbbbqq)
    #149 done
        hist_vec_HHZ_parton_list.append(h_dalpha_rj6_gjs_bbbbqq)
        
        hist_vec_HHZ_parton_list.append(h_mass_comb_gj1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_comb_gj2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_comb_gj3_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_comb_rj1_HHZ_bbbbqq)
    #154 DONE
        hist_vec_HHZ_parton_list.append(h_mass_comb_rj2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_comb_rj3_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj1_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj2_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj3_gjs_bbbbqq)
    #159 done
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj1_Hs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj1_Z_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj2_Hs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj2_Z_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj3_Hs_bbbbqq)
    #164 done
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj3_Z_bbbbqq)
        
        hist_vec_HHZ_parton_list.append(h_mass_comb_rfj_rj1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_comb_rfj_rj2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_comb_rfj_rj3_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rfj_rj1_gjs_bbbbqq)
    #169 done
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rfj_rj2_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rfj_rj3_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rfj_rj1_Hs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rfj_rj1_Z_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rfj_rj2_Hs_bbbbqq)
    #174 done
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rfj_rj2_Z_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rfj_rj3_Hs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_rfj_rj3_Z_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagSum_comb_rfj_rj1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagSum_comb_rfj_rj2_HHZ_bbbbqq)
    #179 done
        hist_vec_HHZ_parton_list.append(h_BTagSum_comb_rfj_rj3_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_comb_rfj_rj1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_comb_rfj_rj2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_comb_rfj_rj3_HHZ_bbbbqq)
        
        hist_vec_HHZ_parton_list.append(h_BTagMax_rfj1_BTag_HHZ_bbbbqq)
    #184 done
        hist_vec_HHZ_parton_list.append(h_BTagMax_rfj2_BTag_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_rfj3_BTag_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_rfj4_BTag_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_rfj5_BTag_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_rfj6_BTag_HHZ_bbbbqq)
    #189 done
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_low_BTag_highMass_Hs_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_low_BTag_highMass_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_low_BTag_highMass_b_Hs_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_low_BTag_highMass_q_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_high_BTag_highMass_Hs_HHZ_bbbbqq)
    #194 done
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_high_BTag_highMass_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_high_BTag_highMass_b_Hs_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_high_BTag_highMass_q_Z_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_low_BTag_highMass_gj_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_rfj_high_BTag_highMass_gj_HHZ_bbbbqq)
    #199 done
        hist_vec_HHZ_parton_list.append(h_mass_comb_rj1_HHZ_bbbbqq_forRMS90)
        hist_vec_HHZ_parton_list.append(h_mass_comb_rj2_HHZ_bbbbqq_forRMS90)
        hist_vec_HHZ_parton_list.append(h_mass_comb_rj3_HHZ_bbbbqq_forRMS90)
        hist_vec_HHZ_parton_list.append(h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90)
        hist_vec_HHZ_parton_list.append(h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90)
    #204 done
        hist_vec_HHZ_parton_list.append(h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90)
        hist_vec_HHZ_parton_list.append(h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90)
        hist_vec_HHZ_parton_list.append(h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90)
        hist_vec_HHZ_parton_list.append(h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90)
        hist_vec_HHZ_parton_list.append(h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq)
    #209 done
        hist_vec_HHZ_parton_list.append(h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_wBTag_rfj_rj1_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_wBTag_rfj_rj2_gjs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_wBTag_rfj_rj3_gjs_bbbbqq)
    #214 done
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_wBTag_rfj_rj1_Hs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_wBTag_rfj_rj1_Z_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_wBTag_rfj_rj2_Hs_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_wBTag_rfj_rj2_Z_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_wBTag_rfj_rj3_Hs_bbbbqq)
    #219 done
        hist_vec_HHZ_parton_list.append(h_dalpha_comb_wBTag_rfj_rj3_Z_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_comb_wBTag_rfj_rj1_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_comb_wBTag_rfj_rj2_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_BTagMax_comb_wBTag_rfj_rj3_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_comb_rj1_matchPart_P_reco_over_P_part_HHZ_bbbbqq)
    #224 done
        hist_vec_HHZ_parton_list.append(h_comb_rj2_matchPart_P_reco_over_P_part_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_comb_rj3_matchPart_P_reco_over_P_part_HHZ_bbbbqq)
    #226 done
        hist_vec_HHZ_parton_list.append(h_comb_rfj_rj1_matchPart_P_reco_over_P_part_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_comb_rfj_rj2_matchPart_P_reco_over_P_part_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_comb_rfj_rj3_matchPart_P_reco_over_P_part_HHZ_bbbbqq)

        hist_vec_HHZ_parton_list.append(h_comb_wBTag_rfj_rj1_matchPart_P_reco_over_P_part_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_comb_wBTag_rfj_rj2_matchPart_P_reco_over_P_part_HHZ_bbbbqq)
    #231 done
        hist_vec_HHZ_parton_list.append(h_comb_wBTag_rfj_rj3_matchPart_P_reco_over_P_part_HHZ_bbbbqq)

        hist_vec_HHZ_parton_list.append(h_H1_comb_rj_match_P_reco_over_P_part_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_H2_comb_rj_match_P_reco_over_P_part_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_Z_comb_rj_match_P_reco_over_P_part_HHZ_bbbbqq)

        hist_vec_HHZ_parton_list.append(h_H1_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq)
    #236 done
        hist_vec_HHZ_parton_list.append(h_H2_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_Z_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq)

        hist_vec_HHZ_parton_list.append(h_H1_comb_wBTag_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_H2_comb_wBTag_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_Z_comb_wBTag_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq)
    #241 done

        hist_vec_HHZ_parton_list.append(h_dalpha_H1_comb_rj_match_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_H2_comb_rj_match_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_Z_comb_rj_match_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_H1_comb_rfj_rj_match_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_H2_comb_rfj_rj_match_HHZ_bbbbqq)
    #246 done
        hist_vec_HHZ_parton_list.append(h_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq)

        hist_vec_HHZ_parton_list.append(h_dalpha_H1_comb_wBTag_rfj_rj_match_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_H2_comb_wBTag_rfj_rj_match_HHZ_bbbbqq)
        hist_vec_HHZ_parton_list.append(h_dalpha_Z_comb_wBTag_rfj_rj_match_HHZ_bbbbqq)
    #250 done

        for hist in hist_vec_HHZ_parton_list:
            hist.Sumw2()

    #print 'low edges of 90', h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.GetBinLowEdge(180),h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.GetBinLowEdge(322)
    #print 'low edges of Z', h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.GetBinLowEdge(112),h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.GetBinLowEdge(254)


    print 'length of     hist list ',len(hist_vec_HHZ_parton_list)
    fill_HHZ_histograms(input_file_,xsec_,_mytree,hist_vec_HHZ_parton_list,lumi,i_ishhzfile)
  
  
    resolution=0.
    resolutionError=0.
    mean=0. 
    meanError=0.
    if i_fillPartonHistos:
        print 'rj1 values 35 window/total def comb', h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.Integral(180,322)/h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj1 values 35 window/total rfj def comb', h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(180,322)/h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj1 values 35 window/total rfj BTag comb', h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(180,322)/h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        resolution, resolutionError, mean, meanError=CalculatePerformance(h_mass_comb_rj1_HHZ_bbbbqq_forRMS90)
        print 'rj1 def comb mean90/err',mean,meanError,'rms90',resolution, resolutionError
        resolution, resolutionError, mean, meanError=CalculatePerformance(h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90)
        print 'rfj_rj1 def comb mean90/err',mean,meanError,'rms90',resolution, resolutionError
        resolution, resolutionError, mean, meanError=CalculatePerformance(h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90)
        print 'rfj_rj1 wBTag comb mean90/err',mean,meanError,'rms90',resolution, resolutionError
        
        print 'rj1 values 35 window/total def comb low', h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.Integral(0,180)/h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj1 values 35 window/total rfj def comb low', h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(0,180)/h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj1 values 35 window/total rfj BTag comb low', h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(0,180)/h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        
        print 'rj1 values 35 window/total def comb high', h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.Integral(322,1501)/h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj1 values 35 window/total rfj def comb high', h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(322,1501)/h_mass_comb_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj1 values 35 window/total rfj BTag comb high', h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(322,1501)/h_mass_comb_wBTag_rfj_rj1_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        
        print 'rj2 values 35 window/total def comb', h_mass_comb_rj2_HHZ_bbbbqq_forRMS90.Integral(180,322)/h_mass_comb_rj2_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj2 values 35 window/total rfj def comb', h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(180,322)/h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj2 values 35 window/total rfj BTag comb', h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(180,322)/h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        resolution, resolutionError, mean, meanError=CalculatePerformance(h_mass_comb_rj2_HHZ_bbbbqq_forRMS90)
        print 'rj2 def comb mean90/err',mean,meanError,'rms90',resolution, resolutionError
        resolution, resolutionError, mean, meanError=CalculatePerformance(h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90)
        print 'rfj_rj2 def comb mean90/err',mean,meanError,'rms90',resolution, resolutionError
        resolution, resolutionError, mean, meanError=CalculatePerformance(h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90)
        print 'rfj_rj2 wBTag comb mean90/err',mean,meanError,'rms90',resolution, resolutionError
        
        print 'rj2 values 35 window/total def comb low', h_mass_comb_rj2_HHZ_bbbbqq_forRMS90.Integral(0,180)/h_mass_comb_rj2_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj2 values 35 window/total rfj def comb low', h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(0,180)/h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj2 values 35 window/total rfj BTag comb low', h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(0,180)/h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        
        print 'rj2 values 35 window/total def comb high', h_mass_comb_rj2_HHZ_bbbbqq_forRMS90.Integral(322,1501)/h_mass_comb_rj2_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj2 values 35 window/total rfj def comb high', h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(322,1501)/h_mass_comb_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj2 values 35 window/total rfj BTag comb high', h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(322,1501)/h_mass_comb_wBTag_rfj_rj2_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        
        print 'rj3 values 35 window/total def comb', h_mass_comb_rj3_HHZ_bbbbqq_forRMS90.Integral(112,254)/h_mass_comb_rj3_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj3 values 35 window/total rfj def comb', h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(112,254)/h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj3 values 35 window/total rfj BTag comb', h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(112,254)/h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        resolution, resolutionError, mean, meanError=CalculatePerformance(h_mass_comb_rj3_HHZ_bbbbqq_forRMS90)
        print 'rj3 def comb mean90/err',mean,meanError,'rms90',resolution, resolutionError
        resolution, resolutionError, mean, meanError=CalculatePerformance(h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90)
        print 'rfj_rj3 def comb mean90/err',mean,meanError,'rms90',resolution, resolutionError
        resolution, resolutionError, mean, meanError=CalculatePerformance(h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90)
        print 'rfj_rj3 wBTag comb mean90/err',mean,meanError,'rms90',resolution, resolutionError
        
        print 'rj3 values 35 window/total def comb low', h_mass_comb_rj3_HHZ_bbbbqq_forRMS90.Integral(0,112)/h_mass_comb_rj3_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj3 values 35 window/total rfj def comb low', h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(0,112)/h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj3 values 35 window/total rfj BTag comb low', h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(0,112)/h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        
        print 'rj3 values 35 window/total def comb high', h_mass_comb_rj3_HHZ_bbbbqq_forRMS90.Integral(254,1501)/h_mass_comb_rj3_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj3 values 35 window/total rfj def comb high', h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(254,1501)/h_mass_comb_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        print 'rj3 values 35 window/total rfj BTag comb high', h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(254,1501)/h_mass_comb_wBTag_rfj_rj3_HHZ_bbbbqq_forRMS90.Integral(0,1501)
        
        
        print 'low edges of Z', h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.GetBinLowEdge(112),h_mass_comb_rj1_HHZ_bbbbqq_forRMS90.GetBinLowEdge(254)

    CLICdpStyle()


    file_histogram.Write()
    file_histogram.Close()
  
    return None

def process_files():

    lumi_=4000.
    print 'start processing of files'

    fillPartonHistos_=True
    ishhzfile_=False
#9,8,5
    cross_section_= 4.18E-02 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_hhqq_14364_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_hhqq_14364_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    ishhzfile_=True
    cross_section_= 6.06e-02
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    fillPartonHistos_=False

    cross_section_= 3.83
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_hzqq_13391_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_hzqq_13391_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    cross_section_= 1269.
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_qq_13399_to_13402_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_qq_13399_to_13402_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    cross_section_= 902.
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_qqqq_13394_to_13397_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_qqqq_13394_to_13397_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    ##here this one has to be run, once I realize the qqqq dataset is finished
    cross_section_= 9.2271753E-03
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_bbcbbc_13094_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_bbcbbc_13094_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    ##here and lower all things should have been run by now
    cross_section_= 9.1731760E-03 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_bbubbu_13095_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_bbubbu_13095_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 1.3757137E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_ddcyyc_13096_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_ddcyyc_13096_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 1.4498909E+01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_dduyyu_13097_polm80_3TeV_wO_CLIC_o3_v14.root" 
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_dduyyu_13097_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 1.2499614E+01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_sscbbc_13098_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_sscbbc_13098_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 1.1651315E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_sscssc_13099_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_sscssc_13099_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 1.2615661E-02 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_ssussu_13123_polm80_3TeV_wO_CLIC_o3_v14.root" 
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_ssussu_13123_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_=  5.4145233E-02 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_ssubbu_13292_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_ssubbu_13292_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_=  1.3394883E+01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_yycbbu_13318_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_yycbbu_13318_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_   

    cross_section_=  2.0054737E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_yycddu_13326_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_yycddu_13326_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_ 

    cross_section_= 2.0248353E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_yycssu_13323_polm80_3TeV_wO_CLIC_o3_v14.root" 
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_yycssu_13323_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_ 

    cross_section_=  1.3330064E+01
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_yyubbc_13320_polm80_3TeV_wO_CLIC_o3_v14.root" 
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_yyubbc_13320_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_=  2.0034170E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_yyuddc_13328_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_yyuddc_13328_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_=  2.0189010E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_yyussc_13325_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_yyussc_13325_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_               
    
    cross_section_=  2.0189010E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_yyussc_13325_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_yyussc_13325_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_                

    cross_section_=  3.6427E-02 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_bbxxxx_14601_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_bbxxxx_14601_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  7.8948E-03 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_bdxxxx_14603_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_bdxxxx_14603_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  7.8901E-03
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_bsxxxx_14604_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_bsxxxx_14604_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=   7.9020E-03
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_dbxxxx_14605_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_dbxxxx_14605_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_= 2.3698E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_ddxxxx_14481_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_ddxxxx_14481_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  4.3881E+00
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_dsxxxx_14606_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_dsxxxx_14606_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=   8.0580E-03 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_sbxxxx_14485_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_sbxxxx_14485_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_= 4.4150E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_sdxxxx_14487_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_sdxxxx_14487_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  2.2316E+00
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polm80/HHZStudy_ee_ssxxxx_14609_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_ssxxxx_14609_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    lumi_=1000.
    fillPartonHistos_=True
    ishhzfile_=False
#9,8,5
    cross_section_=   2.898E-02 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_hhqq_14365_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_hhqq_14365_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    ishhzfile_=True
    cross_section_=  4.23E-02 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_hhz_14344_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_hhz_14344_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    fillPartonHistos_=False

    cross_section_=  2.67
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_hzqq_13392_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_hzqq_13392_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    cross_section_= 786.
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_qq_13398_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_qq_13398_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    cross_section_= 120.
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_qqqq_13393_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_qqqq_13393_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    cross_section_= 2.9986901E-03 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_bbcbbc_13071_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_bbcbbc_13071_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_)
    print 'finished file', final_histo_name_

    cross_section_= 2.9825397E-03
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_bbubbu_13072_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_bbubbu_13072_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 1.7824610E-01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_ddcyyc_13073_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_ddcyyc_13073_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 5.0109474E+00
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_dduyyu_13074_polp80_3TeV_wO_CLIC_o3_v14.root" 
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_dduyyu_13074_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 4.8938333E+00
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_sscbbc_13075_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_sscbbc_13075_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 1.3677677E-01
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_sscssc_13076_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_sscssc_13076_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 3.3776171E-03
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_ssussu_13077_polp80_3TeV_wO_CLIC_o3_v14.root" 
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_ssussu_13077_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_=  2.3216638E-02
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_ssubbu_13293_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_ssubbu_13293_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 5.2101109E+00
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_yycbbu_13322_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_yycbbu_13322_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_   

    cross_section_= 4.0984879E-01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_yycddu_13319_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_yycddu_13319_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_ 

    cross_section_= 4.1853929E-01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_yycssu_13327_polp80_3TeV_wO_CLIC_o3_v14.root" 
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_yycssu_13327_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_ 

    cross_section_= 5.2070149E+00 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_yyubbc_13324_polp80_3TeV_wO_CLIC_o3_v14.root" 
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_yyubbc_13324_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 4.1203686E-01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_yyuddc_13321_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_yyuddc_13321_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_

    cross_section_= 4.2245034E-01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_yyussc_13329_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/ntuple_HHZ_ee_yyussc_13329_polm80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_               

    cross_section_= 1.1996E-02 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_bbxxxx_14602_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_bbxxxx_14602_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  9.1650E-04
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_bdxxxx_14476_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_bdxxxx_14476_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=   9.0803E-04
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_bsxxxx_14478_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_bsxxxx_14478_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  9.1236E-04 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_dbxxxx_14480_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_dbxxxx_14480_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  2.7583E-01
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_ddxxxx_14482_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_ddxxxx_14482_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  5.0859E-01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_dsxxxx_14607_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_dsxxxx_14607_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=   9.1615E-04 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_sbxxxx_14608_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_sbxxxx_14608_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  5.1915E-01  
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_sdxxxx_14488_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_sdxxxx_14488_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  

    cross_section_=  2.6630E-01 
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC7_NJet6_finalAnalysis/polp80/HHZStudy_ee_ssxxxx_14610_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/ntuple_HHZ_ee_ssxxxx_14610_polp80_3TeV_wO_CLIC_o3_v14.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,fillPartonHistos_,ishhzfile_) 
    print 'finished file', final_histo_name_  


    return None

process_files()



