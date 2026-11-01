from django.http import HttpResponse, HttpResponseRedirect
from django.db.models import Q
from django.shortcuts import render_to_response
from models import all, bestMX
import settings
#import settings
from common import *
import re

# Create your views here.
def search_form(request):
    return render_to_response('search_form.html')

def example(request, MXorSeq):
  errors = []
  if "seq" in MXorSeq:
    q = "TGACTC"
    return render_to_response('form.html', {'seq': q, 'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
  if "mx" in MXorSeq:
    r = "A | 0.515 0.030 0.064 0.687 0.044 0.117 0.382\nC | 0.087 0.028 0.885 0.056 0.561 0.182 0.178\nG | 0.317 0.045 0.024 0.190 0.024 0.164 0.138\nT | 0.081 0.897 0.027 0.067 0.371 0.537 0.303"
    return render_to_response('form.html', {'mx': r, 'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
  errors.append("Error in Loading Example Data")
  return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

def matrix(request):
  string = request.META['QUERY_STRING']
  vec = string.split('+')
  db = vec[0]
  dbmod = vec[1]
  tf = vec[-1].split('=')
  name = tf[0].split('_')
  mxs = all.objects.filter( Q(name=tf[0]) & Q(source__icontains=db) & Q(source__icontains=dbmod))
  return render_to_response('matrix.html', {'mxs': mxs, 'mx_type': tf[1], 'tf_name': name[0], 'STATIC_URL' : settings.STATIC_URL, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

def details(request, TFname):
  errors = []
  q = TFname
  q = re.sub(r'\s','',q)
  r = q + "_"

  bests=bestMX.objects.filter(name=q)
  alls = all.objects.filter( Q(name=q) | Q(name__icontains=r))
  mxs=[]
  for Bmx in bests:
    for Amx in alls:
      if Bmx.source == Amx.source and Bmx.file == Amx.file: continue
      mxs.append(Amx)

  if len(bests)==0:
    errors.append('Unable to find TF in database')
    return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

  sorted_mxs = sorted(mxs, key=lambda x:x.pvalue)
  for mx in sorted_mxs:
    if mx.pvalue:
      if mx.pvalue < 0.0001: mx.pvalue = "%.2e" % mx.pvalue
      else: mx.pvalue = "%1.4f" % mx.pvalue

  for Bmx in bests:
    if len(Bmx.chip)>200:
      genes = Bmx.chip.split(", ")
      Bmx.chip = ", ".join(genes[:18]) + ' <a href="/ScerTF/CHIP/'+Bmx.name+'">+' + str(len(genes)-18) +' more</a>'
    if len(Bmx.knockout)>200:
      genes = Bmx.knockout.split(", ")
      Bmx.knockout = ", ".join(genes[:17]) + ' <a href="/ScerTF/KO/'+Bmx.name+'">+' + str(len(genes)-17) +' more</a>'
    if Bmx.pvalue:
      if Bmx.pvalue < 0.0001: Bmx.pvalue = "%.2e" % Bmx.pvalue
      else: Bmx.pvalue = "%1.4f" % Bmx.pvalue
  return render_to_response('details.html', {'bests': bests, 'mxs': sorted_mxs, 'query': q, 'STATIC_URL' : settings.STATIC_URL, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

def systematic(request, TFname):
  errors = []
  q = "Y" + TFname
  q = re.sub(r'\s','',q)

  bests=bestMX.objects.filter(alias__icontains=q)
  if len(bests) >= 1:
    for Bmx in bests:
	q = Bmx.name

  r = q + "_"
  bests=bestMX.objects.filter(name=q)
  alls = all.objects.filter( Q(name=q) | Q(name__icontains=r))
  mxs=[]
  for Bmx in bests:
    for Amx in alls:
      if Bmx.source == Amx.source and Bmx.file == Amx.file: continue
      mxs.append(Amx)

  if len(bests)==0:
    errors.append('Unable to find TF in database')
    return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

  sorted_mxs = sorted(mxs, key=lambda x:x.pvalue)
  for mx in sorted_mxs:
    if mx.pvalue:
      if mx.pvalue < 0.0001: mx.pvalue = "%.2e" % mx.pvalue
      else: mx.pvalue = "%1.4f" % mx.pvalue

  for Bmx in bests:
    if len(Bmx.chip)>200:
      genes = Bmx.chip.split(", ")
      Bmx.chip = ", ".join(genes[:18]) + ' <a href="/ScerTF/CHIP/'+Bmx.name+'">+' + str(len(genes)-18) +' more</a>'
    if len(Bmx.knockout)>200:
      genes = Bmx.knockout.split(", ")
      Bmx.knockout = ", ".join(genes[:17]) + ' <a href="/ScerTF/KO/'+Bmx.name+'">+' + str(len(genes)-17) +' more</a>'
    if Bmx.pvalue:
      if Bmx.pvalue < 0.0001: Bmx.pvalue = "%.2e" % Bmx.pvalue
      else: Bmx.pvalue = "%1.4f" % Bmx.pvalue
  return render_to_response('details.html', {'bests': bests, 'mxs': sorted_mxs, 'query': q, 'STATIC_URL' : settings.STATIC_URL, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

def search(request):
    nt = ('A', 'C', 'G', 'T')
    alpha = 0.05
    errors = []
    if 'clear' in request.GET:
      return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
    if 'q' in request.GET:
        q = request.GET['q']
        q = re.sub(r'\s','',q)
        if q == 'sailboat': 
          return HttpResponseRedirect("http://i.imgur.com/1vFrP.gif")
        if q == 'vitamins': 
          return HttpResponseRedirect("http://i.imgur.com/XJyGY.gif")
        if not q: 
            errors.append('Please enter a search term.')
        elif len(q) > 20:
            errors.append('Error: input seqeuence over 20 characters. If entering a TF name or DNA sequence, please use at most 20 characters.')
        if q:
            mxs = bestMX.objects.filter(name=q)
            if mxs:
              return HttpResponseRedirect("/ScerTF/details/"+q+"/")
            if(len(q)<2):
              errors.append('Error: Query Too Short. Use at least 2 characters.')
              return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
            if not mxs: mxs = bestMX.objects.filter(name__icontains=q)
            if not mxs: mxs = bestMX.objects.filter(alias__icontains=q)
            if not mxs:
                if(re.search('[^ACGTacgt]', q)):
                  errors.append('Error: No TF gene names matched your search, and your query contained non-A,C,G,T letters')
                  return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
                if len(q) < 4:
                  errors.append('Error: No TF gene names matched your search, so I assumed it was a nucleotide sequence.  Nucleotide sequences must be at least 4 long.')
                  return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
                r = SeqToMX(q)
                mxs = bestMX.objects.filter(name__icontains='')
                scores = []
                d = {}
                FullSite = {}
                for mx in mxs:
                  Fscore = align(r,mx.PFM)
                  rev=reverse_MX(r)
                  Rscore = align(rev,mx.PFM)
                  if Fscore >= Rscore:
                    d[mx.name] = seqZtest(q,mx.PFM)
                    scores.append(d[mx.name]);
                    FullSite[mx.name] = Fscore[4]-Fscore[3]+1
                  else:
                    d[mx.name] = seqZtest(reverse_complement(q),mx.PFM)
                    scores.append(d[mx.name]);
                    FullSite[mx.name] = Rscore[4]-Rscore[3]+1
                hits = sortPvals(scores,d)
                mxs = []
                for TF in hits:
                  tmp=bestMX.objects.filter(name=TF)
                  for mx in tmp:
                    mx.pvalue = hits[TF]
                    if len(mx.consensus) > FullSite[mx.name]: mx.note = "Partial Overlap"
                    else: mx.note = "NA"
                    mxs.append(mx)
                sorted_mxs = sorted(mxs, key=lambda x:x.pvalue)
                for mx in sorted_mxs:
                  if mx.pvalue < 0.0001: mx.pvalue = "%.2e" % mx.pvalue
                  else: mx.pvalue = "%1.4f" % mx.pvalue
               #assert false
                return render_to_response('search_results.html', 
                      {'mxs': sorted_mxs, 'query': q, 'STATIC_URL' : settings.STATIC_URL, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
            for mx in mxs: mx.pvalue="NA"
            return render_to_response('search_results.html', 
                  {'mxs': mxs, 'query': q, 'STATIC_URL' : settings.STATIC_URL, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
    if 'r' in request.GET:
        r = request.GET['r']
        if r: 
          #Validate Input Matrix.  Is it actually a matrix?  Is it a PFM, PWM, or PCM?
          r=r.upper()
          r = re.sub(r'[^\d\n.ACGT\s\-\|]+','',r)
          if not re.search('A\s\|', r):  r = re.sub(r'A','A |',r)
          if not re.search('C\s\|', r):  r = re.sub(r'C','C |',r)
          if not re.search('G\s\|', r):  r = re.sub(r'G','G |',r)
          if not re.search('T\s\|', r):  r = re.sub(r'T','T |',r)
          if (r.count("A") > 1) or (r.count("C") > 1) or (r.count("G") > 1) or (r.count("T") > 1):
           errors.append('Matrix contains more than one entry for either A, C, G, or T')
           return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
          if len(r.split('\n')) < len(nt): 
           errors.append("Matrix Input Error:  Not enough rows. Nucleotides are organized as rows. Each column is a different position in the matrix")
           return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
          m = [line.split('|') for line in r.split('\n')]
          col_sum=0
          for i in range(len(nt)):
           for j in range(len(m)):
            if nt[i] == m[j][0].rstrip():
              row = map(float,m[j][1].split())
              col_sum += row[0]
          if(re.search('-',r)): r = pwm2pfm(r)
          elif(col_sum > 1.1): r = pcm2pfm(r)

          #Search database for similar matrices
          seq=consensus(r)
          mxs = bestMX.objects.filter(name__icontains='')
          scores = []
          d = {}
          FullSite = {}
          for mx in mxs:
            Fscore = align(r,mx.PFM)
            rev=reverse_MX(r)
            Rscore = align(rev,mx.PFM)
            if Fscore >= Rscore:
              d[mx.name] = mxZtest(r,mx.PFM)
              scores.append(d[mx.name]);
              FullSite[mx.name] = Fscore[4]-Fscore[3]+1
              #scores.append(Fscore[0]);
              #d[mx.name] = Fscore[0]
            else:
              d[mx.name] = mxZtest(rev,mx.PFM)
              scores.append(d[mx.name]);
              FullSite[mx.name] = Rscore[4]-Rscore[3]+1
              #scores.append(Rscore[0]);
              #d[mx.name] = Rscore[0]
          #hits = Ztest(scores,d)
          hits = sortPvals(scores,d)
          mxs = []
          for TF in hits:
            tmp=bestMX.objects.filter(name=TF)
            for mx in tmp:
              mx.pvalue = hits[TF]
              if len(mx.consensus) > FullSite[mx.name]: mx.note = "Partial Overlap"
              else: mx.note = "NA"
              mxs.append(mx)
          sorted_mxs = sorted(mxs, key=lambda x:x.pvalue)
          if sorted_mxs[0].pvalue > alpha : 
            sorted_mxs[0].pvalue = sorted_mxs[0].pvalue / len(scores)
            sorted_mxs[0].note = "No Significant Hits"
          #assert false
          for mx in sorted_mxs:
            if mx.pvalue < 0.0001: mx.pvalue = "%.2e" % mx.pvalue
            else: mx.pvalue = "%1.4f" % mx.pvalue
          #assert false
          return render_to_response('search_results.html', 
                {'mxs': sorted_mxs, 'query': seq, 'STATIC_URL' : settings.STATIC_URL, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
    return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

def download(request):
  return render_to_response('download.html',{'TEMPLATE_URL' : settings.TEMPLATE_URL, 'STATIC_URL' : settings.STATIC_URL})

def about(request):
  return render_to_response('about.html',{'TEMPLATE_URL' : settings.TEMPLATE_URL})

def references(request):
  return render_to_response('data_sources.html',{'TEMPLATE_URL' : settings.TEMPLATE_URL})

def browse(request):
  mxs = bestMX.objects.filter(name__icontains='').order_by('name')
  return render_to_response('browse.html',{'mxs': mxs, 'TEMPLATE_URL' : settings.TEMPLATE_URL, 'STATIC_URL' : settings.STATIC_URL})


def genes(request, TFname):
  errors = []
  q = TFname
  q = re.sub(r'\s','',q)
  r = q + "_"

  bests=bestMX.objects.filter(name=q)
  if len(bests)==0:
    errors.append('Unable to find TF in database')
    return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

  for Bmx in bests:
    both = set(Bmx.chip.split(", ")).intersection(set(Bmx.knockout.split(", ")))
    genes = ", ".join(both)
  return render_to_response('genes.html', {'bests': bests, 'query': q, 'intersect': genes, 'STATIC_URL' : settings.STATIC_URL, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

def pairs(request, TFname):
  errors = []
  q = TFname
  q = re.sub(r'\s','',q)
  r = q + "_"

  bests=bestMX.objects.filter(name=q)
  if len(bests)==0:
    errors.append('Unable to find TF in database')
    return render_to_response('form.html', {'errors': errors, 'TEMPLATE_URL' : settings.TEMPLATE_URL})

  for Bmx in bests:
    tfs = map(lambda x: Bmx.name+"x"+x,Bmx.pairs.split(", "))
  return render_to_response('pairs.html', {'bests': bests, 'pairs': tfs, 'query': q, 'STATIC_URL' : settings.STATIC_URL, 'TEMPLATE_URL' : settings.TEMPLATE_URL})
