from django.db import models

# Create your models here.
class combinedMX(models.Model):
    name = models.CharField(max_length=40)
    consensus = models.CharField(max_length=40)
    PFM = models.CharField(max_length=2000)
    logofile = models.CharField(max_length=40)
    #pvalue = models.FloatField()

    def __unicode__(self):
        return self.name


class all(models.Model):
    name = models.CharField(max_length=40)
    file = models.CharField(max_length=40)
    source = models.CharField(max_length=40)
    ref = models.CharField(max_length=100)
    consensus = models.CharField(max_length=40)
    PFM = models.CharField(max_length=2000)
    PCM = models.CharField(max_length=2000)
    PWM = models.CharField(max_length=2000)
    logofile = models.CharField(max_length=40)
    pvalue = models.FloatField(blank=True, null=True)
    cutoff = models.FloatField()
    ic = models.FloatField()
    max = models.FloatField()
    def __unicode__(self):
        return self.name

class bestMX(models.Model):
    name = models.CharField(max_length=40)
    file = models.CharField(max_length=40)
    source = models.CharField(max_length=40)
    ref = models.CharField(max_length=100)
    alias = models.CharField(max_length=250)
    consensus = models.CharField(max_length=40)
    PFM = models.CharField(max_length=2000)
    logofile = models.CharField(max_length=40)
    description = models.TextField()
    pvalue = models.FloatField(blank=True, null=True)
    cutoff = models.FloatField()
    ic = models.FloatField()
    max = models.FloatField()
    chip = models.TextField()
    knockout = models.TextField()
    pairs = models.TextField()
    note = models.TextField()

    def __unicode__(self):
        return self.name
