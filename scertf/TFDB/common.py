#Common python modules that will be useful in the project
#from numpy import zeros
from math import sqrt,exp,log
import re

def LogOdds(x, y): return (log(float(x)/float(y)))

def reverse_complement(s):
  comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
  #s = re.sub('[^ACGTacgt]', '', s)
  if(re.search('[^ACGTacgt]', s)): return 'Error: Input contains non-A,C,G,T letters'
  letters = list(s.upper())
  letters = [comp[base] for base in letters]
  letters.reverse()
  return ''.join(letters)

def reverse_MX(a):
  lines_a = a.split('\n')
  m_a = [line.split('|') for line in lines_a]
  m = ["A |", "C |", "G |", "T |"]
  nt = ('T', 'G', 'C', 'A')
  for i in range(len(nt)):
    for j in range(len(m_a)):
      if nt[i] == m_a[j][0].rstrip():
        m_a[j][1] = re.sub('[^\d. ]', '',  m_a[j][1])
        temp = map(float,m_a[j][1].split())
        m[i] = m[i] + " ".join(map(str,temp[::-1]))
  matrix = m[0] + " \n" + m[1] + " \n" + m[2] + " \n" + m[3] + " \n"
  return matrix


def consensus(a):
  lines_a = a.split('\n')
  m_a = [line.split('|') for line in lines_a]

  nt = ('A', 'C', 'G', 'T')
  m1 = []
  for i in nt:
    for j in range(len(m_a)):
      if i == m_a[j][0].rstrip():
        m_a[j][1] = re.sub('[^\d. ]', '',  m_a[j][1])
        m1.append(map(float,m_a[j][1].split()))
  if len(m1) < len(nt):
    return -1

  seq = ""
  for i in range(0,len(m1[0])):
    col = map(float,[row[i] for row in m1])
    val =  col.index(max(col))
    seq = seq + nt[val]
  return seq


def SeqToMX(s):
  s = s.upper()
  s = re.sub('[^ACGT]', '', s)
  nt = ['A','C','G','T']
  m = ["A |", "C |", "G |", "T |"]
  for i in range(0,len(s)):
    val = nt.index(s[i])
    for j in range(0,len(m)):
      if j == val:
        m[j] = m[j] + " 0.85"
      else:
        m[j] = m[j] + " 0.05"
  matrix = m[0] + " \n" + m[1] + " \n" + m[2] + " \n" + m[3] + " \n"
  return matrix
    
def pwm2pfm(input):
  lines_a = input.split('\n')
  m_a = [line.split('|') for line in lines_a]

  nt = ('A', 'C', 'G', 'T')
  m1 = []
  for i in nt:
    for j in range(len(m_a)):
      if i == m_a[j][0].rstrip():
        m_a[j][1] = re.sub('[^\d. \-]', '',  m_a[j][1])
        m1.append(map(float,m_a[j][1].split()))
  if len(m1) < len(nt):
    return -1

  #Compute Background frequencies for ACGT from PWM
  bgATtot = 0
  failures = 0
  for j in range(0,len(m1[0])):
   lAT = exp(m1[0][j]) + exp(m1[3][j])
   lCG = exp(m1[1][j]) + exp(m1[2][j])
   bgAT = (1.0 - 0.5*lCG)/(lAT-lCG)
   if bgAT > 1 or bgAT < 0:
     failures+=1
     continue
   bgATtot += bgAT
  m1_len = len(m1[0]) - failures
  bgAT = bgATtot / m1_len
  bgCG = 0.5 - bgAT
  
  #Using background frequency calculated above, calculate PFM
  for j in range(0, len(m1[0])):
    m1[0][j] = exp(m1[0][j])*bgAT
    m1[1][j] = exp(m1[1][j])*bgCG
    m1[2][j] = exp(m1[2][j])*bgCG
    m1[3][j] = exp(m1[3][j])*bgAT

  #convert matrix back to character string
  matrix = "A | " + " ".join(map(str,(m1[0]))) + " \n" + "C | " + " ".join(map(str,(m1[1]))) + " \n" + "G | " + " ".join(map(str,(m1[2]))) + " \n" + "T | " + " ".join(map(str,(m1[3]))) + " \n"
  return matrix

def pcm2pfm(input):
  lines_a = input.split('\n')
  m_a = [line.split('|') for line in lines_a]

  nt = ('A', 'C', 'G', 'T')
  m1 = []
  for i in nt:
    for j in range(len(m_a)):
      if i == m_a[j][0].rstrip():
        m_a[j][1] = re.sub('[^\d. ]', '',  m_a[j][1])
        m1.append(map(float,m_a[j][1].split()))
  if len(m1) < len(nt):
    return -1

  for j in range(0,len(m1[0])):
    col = sum(map(float,[row[j] for row in m1]))
    for i in range(0,len(m1)):
      m1[i][j] = m1[i][j] / col
      if(m1[i][j] <0.001): m1[i][j] = 0.001

  #convert matrix back to character string
  matrix = "A | " + " ".join(map(str,(m1[0]))) + " \n" + "C | " + " ".join(map(str,(m1[1]))) + " \n" + "G | " + " ".join(map(str,(m1[2]))) + " \n" + "T | " + " ".join(map(str,(m1[3]))) + " \n"
  return matrix

def pearsonr(x,y):
  n=len(x)
  vals=range(n)
     
  #regular sums
  sumx=sum([float(x[i]) for i in vals])
  sumy=sum([float(y[i]) for i in vals])
     
  #sum of the squares
  sumxSq=sum([x[i]**2.0 for i in vals])
  sumySq=sum([y[i]**2.0 for i in vals])
   
  #sum of the products
  pSum=sum([x[i]*y[i] for i in vals])
   
  #do pearson score
  num=pSum-(sumx*sumy/n)
  denX=(sumxSq-pow(sumx,2)/n) 
  denY=(sumySq-pow(sumy,2)/n)
  den=(denX*denY)
  #print n, " ", vals, " ", sumxSq, " ", sumySq, " ", pSum, " ", num, " ", den, " ", sumx, " ", sumy, " ", denY, " ", pow(sumy,2)
  if den <= 0:
    return 1
  den=sqrt(den)
  r=num/den

  #print den, " ", r

  if r > 1:
    return 1
  elif r < -1:
    return -1
  else:
    return r

def allr(x,y):
  bg = (0.32, 0.18, 0.18, 0.32)
  nt = ('A', 'C', 'G', 'T')
  if(len(x) != len(y)):
    return -10
  sum = 0
  for i in range(0,len(x)):
    sum += (x[i]*(y[i]/bg[i])) + (y[i]*(x[i]/bg[i]))
  sum = sum/2
  return sum

def align(a,b):
  a = re.sub('\t', ' ',a);
  b = re.sub('\t', ' ',b);
  lines_a = a.split('\n')
  lines_b = b.split('\n')
  m_a = [line.split('|') for line in lines_a]
  m_b = [line.split('|') for line in lines_b]

  nt = ('A', 'C', 'G', 'T')

  m1 = []
  for i in nt:
    for j in range(len(m_a)):
      if i == m_a[j][0].rstrip():
        m_a[j][1] = re.sub('[^\d. ]', '',  m_a[j][1])
        m1.append(map(float,m_a[j][1].split()))
  if len(m1) < len(nt):
    return -1

  m2 = []
  for i in nt:
    for j in range(len(m_b)):
      if i == m_b[j][0].rstrip():
        m_b[j][1] = re.sub('[^\d. ]', '',  m_b[j][1])
        m2.append(map(float,m_b[j][1].split()))
  if len(m2) < len(nt):
    return -1

  #Local Smith waterman alignment requiring a minimum alignment length of 4
  #PCC between columns is the scoring metric (PCC ranges -1,1)

  m = len(m1[0])
  n = len(m2[0])
  score = [[0]*(n+1) for x in xrange(m+1)]
  pointer = [[0]*(n+1) for x in xrange(m+1)]
  for i in range(0,m+1): 
    score[i][0] = 0    #the DP table
    pointer[i][0] = 0  #the Traceback table
  for j in range(0,n+1):
    score[0][j] = 0    #the DP table
    pointer[0][j] = 0  #the Traceback table

  #calculate DP table and mark pointers
  for i in range(1,m+1):
    col1 = [row[i-1] for row in m1]
    for j in range(1,n+1):
      col2 = [row[j-1] for row in m2]
      score_diagonal=score[i-1][j-1]+allr(col1,col2)
      score[i][j]=score_diagonal
      if score[i][j]==score_diagonal:
        pointer[i][j]=3 #3 means trace diagonal
  #END of DP table

  #Find max 
  min_length=4
  max_score = -10000 
  max_i=0
  max_j=0
  for i in range(min_length,m+1):
    if score[i][n] > max_score: 
      max_score = score[i][n]
      max_i = i
      max_j = j
  for j in range(min_length,n+1):
    if score[m][j] > max_score: 
      max_score = score[m][j]
      max_i = i
      max_j = j

  #traceback
  mx1_end_align = max_i
  mx2_end_align = max_j
  i = max_i
  j = max_j
  while pointer[i][j]!=0:
    if pointer[i][j]==3:
      i=i-1;
      j=j-1;
  #END of traceback
  mx1_start_align=i+1
  mx2_start_align=j+1
  return (max_score, mx1_start_align, mx1_end_align, mx2_start_align, mx2_end_align)

def Ztest(scores,d):
  hits = {}
  min = 2
  min_name=""
  minDict = {}
  N = len(scores)
  alpha = 0.01
  mean = sum(scores)/N
  SD = sqrt( sum( ((x-mean)**2) for x in scores) / N)
  SE = SD #/ 2 #sqrt(N)
  for name, score in d.iteritems():
    Z = (score-mean)/SE
    pval = (1-pnorm(Z)) #*N
    if pval <= alpha:
      hits[name] = pval
    if pval < min: 
      min = pval
      min_name = name
  if len(hits) == 0:
    minDict[min_name]=min
    return minDict
  return hits

def sortPvals(scores,d):
  hits = {}
  min = 999999999999999
  min_name=""
  minDict = {}
  N = len(scores)
  alpha = 0.01
  for name, pval in d.iteritems():
    pval = pval * N
    if pval <= alpha:
      hits[name] = pval
    if pval < min: 
      min = pval
      min_name = name
  if len(hits) == 0:
    minDict[min_name]=min
    return minDict
  return hits


def seqZtest(s,b):
  s = s.upper()
  nt = ('A', 'C', 'G', 'T')
  bg = (0.32, 0.18, 0.18, 0.32)
  a = SeqToMX(s)
  alignment = align(a,b)

  #Convert the aligned region of the database matrix (b) to a PWM
  start = alignment[3]-1
  end = alignment[4]
  b = re.sub('\t', ' ',b);
  lines_b = b.split('\n')
  m_b = [line.split('|') for line in lines_b]
  db_mx = []
  for i in nt:
    for j in range(len(m_b)):
      if i == m_b[j][0].rstrip():
        m_b[j][1] = re.sub('[^\d. ]', '',  m_b[j][1])
        db_mx.append(map(LogOdds, m_b[j][1].split()[start:end], [bg[list(nt).index(i)]]*(end-start)))
  if len(db_mx) < len(nt):
    return -1

  #Determine mean and standard deviation of scores produced by PWM aligned region
  db_stats = mx_dist(db_mx)

  #Determine score of sequence submitted to subroutine
  start = alignment[1]-1
  end = alignment[2]
  score = 0
  pos = 0
  for i in range(start,end):
    score += db_mx[list(nt).index(s[i])][pos]
    pos+=1
	
  #Given the distribution defined by mean=db_stats[0] and variance=db_stats[1], 
  # what is p-value of getting pwm-score: score
  mean = db_stats[0]
  SD = sqrt(db_stats[1])
  Z = (score-mean)/SD
  pval = (1-pnorm(Z))
  return pval
  
def mx_dist(m):
  mx_mean = 0
  mx_var = 0
  for j in range(0,len(m[0])):
    col_sum = sum(map(float,[row[j] for row in m]))
    col_mean = col_sum / len(m)
    col_var = sum((row[j]-col_mean)**2 for row in m) / len(m)
    mx_mean += col_mean
    mx_var += col_var
  return (mx_mean, mx_var)
 
 

def mxZtest(a,b):
  nt = ('A', 'C', 'G', 'T')
  bg = (0.32, 0.18, 0.18, 0.32)
  alignment = align(a,b)
  
  #Extract the aligned region of the query matrix (PFM)
  start = alignment[1]-1
  end = alignment[2]
  a = re.sub('\t', ' ',a);
  lines_a = a.split('\n')
  m_a = [line.split('|') for line in lines_a]
  query_mx = []
  for i in nt:
    for j in range(len(m_a)):
      if i == m_a[j][0].rstrip():
        m_a[j][1] = re.sub('[^\d. ]', '',  m_a[j][1])
        query_mx.append(map(float,m_a[j][1].split())[start:end])
  if len(query_mx) < len(nt):
    return -1
  
  
  #Convert the aligned region of the database matrix (b) to a PWM
  start = alignment[3]-1
  end = alignment[4]
  b = re.sub('\t', ' ',b);
  lines_b = b.split('\n')
  m_b = [line.split('|') for line in lines_b]
  db_mx = []
  for i in nt:
    for j in range(len(m_b)):
      if i == m_b[j][0].rstrip():
        m_b[j][1] = re.sub('[^\d. ]', '',  m_b[j][1])
        db_mx.append(map(LogOdds, m_b[j][1].split()[start:end], [bg[list(nt).index(i)]]*(end-start)))
  if len(db_mx) < len(nt):
    return -1

  #Determine mean and standard deviation of scores produced by PWM aligned region
  db_stats = mx_dist(db_mx)
  
  #Multiply the user-input PFM by the database PWM to obtain the score of the 
  # average sequence produced by the user-input PFM
  score = 0
  for i in range(len(query_mx[0])):
    for j in range(len(query_mx)):
      score += query_mx[j][i] * db_mx[j][i]

  #Given the distribution defined by mean=db_stats[0] and variance=db_stats[1], 
  # what is p-value of getting pwm-score: score
  mean = db_stats[0]
  SD = sqrt(db_stats[1])
  Z = (score-mean)/SD
  pval = (1-pnorm(Z))
  return pval

def phi(x):
  # constants
  a1 =  0.254829592
  a2 = -0.284496736
  a3 =  1.421413741
  a4 = -1.453152027
  a5 =  1.061405429
  p  =  0.3275911

  # Save the sign of x
  sign = 1
  if x < 0:
      sign = -1
  x = abs(x)/sqrt(2.0)

  # A&S formula 7.1.26
  t = 1.0/(1.0 + p*x)
  y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
  return 0.5*(1.0 + sign*y)

def pnorm(a, mu = 0, sigma = 1):
  return phi((a-mu)/sigma)

def SW(a,b):
  lines_a = a.split('\n')
  lines_b = b.split('\n')
  m_a = [line.split('|') for line in lines_a]
  m_b = [line.split('|') for line in lines_b]

  nt = ('A', 'C', 'G', 'T')

  m1 = []
  for i in nt:
    for j in range(len(m_a)):
      if i == m_a[j][0].rstrip():
        m_a[j][1] = re.sub('[^\d. ]', '',  m_a[j][1])
        m1.append(map(float,m_a[j][1].split()))
  if len(m1) < len(nt):
    return -1

  m2 = []
  for i in nt:
    for j in range(len(m_b)):
      if i == m_b[j][0].rstrip():
        m_b[j][1] = re.sub('[^\d. ]', '',  m_b[j][1])
        m2.append(map(float,m_b[j][1].split()))
  if len(m2) < len(nt):
    return -1

  #Local Smith waterman alignment requiring a minimum alignment length of 4
  #PCC between columns is the scoring metric (PCC ranges -1,1)
  #Forward Pass:
  m = len(m1[0])
  n = len(m2[0])
  score=zeros((m+1,n+1))      #the DP table
  pointer=zeros((m+1,n+1))    #the traceback path
  length = zeros((m+1,n+1))   #the lengths of each path
  
  max_score = 0 
  max_i = max_j = 0
  min_length= 4 
  #calculate DP table and mark pointers
  for i in range(1,m+1):
    col1 = [row[i-1] for row in m1]
    for j in range(1,n+1):
      col2 = [row[j-1] for row in m2]
      score_diagonal=score[i-1][j-1]+allr(col1,col2)
      score[i][j]=max(0,score_diagonal)
      if score[i][j]==0:
        pointer[i][j]=0 #0 means end of the path
        length[i][j]=0
      if score[i][j]==score_diagonal:
        pointer[i][j]=3 #3 means trace diagonal
        length[i][j]=length[i-1][j-1]+1
      if score[i][j]>=max_score and length[i][j]>=min_length:
        max_i=i
        max_j=j
        max_score=score[i][j]
  #END of DP table

  #traceback
  mx1_end_align = max_i
  mx2_end_align = max_j
  i = max_i
  j = max_j
  while pointer[i][j]!=0:
    if pointer[i][j]==3:
      i=i-1;
      j=j-1;
  #END of traceback
  mx1_start_align=i+1
  mx2_start_align=j+1

  return (max_score, mx1_start_align, mx1_end_align, mx2_start_align, mx2_end_align)

 