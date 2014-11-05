## copied verbatim from ngsdb.ngsdbview
# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#     * Rearrange models' order
#     * Make sure each model has one field with primary_key=True
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin.py sqlcustom [appname]'
# into your database.

from django.db import models
from django.contrib.auth.models import User
from ngsdbview.validators import *
#import fields import *

EXPERIMENT_TYPE_CHOICES = (
    ('SL', 'SL'),
    ('RNAseq', 'RNAseq'),
    ('RiboProf', 'RiboProf'),
    ('DNAseq', 'DNAseq'),
    ('mixed', 'mixed')
)


class Author(models.Model):
    author_id = models.AutoField(primary_key=True)
    firstname = models.CharField(max_length=45)
    lastname = models.CharField(max_length=45)
    designation = models.CharField(unique=True, max_length=5)
    email = models.EmailField(unique=True)
    def __unicode__(self):
        return str(self.designation)

class Librarytype(models.Model):
    librarytype_id = models.AutoField(primary_key=True)
    type = models.CharField(unique=True, max_length=25)
    notes = models.TextField(max_length=400, blank=True)
    def __unicode__(self):
        return str(self.type)

class Protocol(models.Model):
    protocol_id = models.AutoField(primary_key=True)
    protocol_name = models.CharField(unique=True, max_length=50)
    protocol_file = models.FileField(upload_to="protocols")
    protocol_link = models.URLField(max_length=1000, blank=True)
    notes = models.CharField(max_length=400, default=None)

    class Meta:
        ordering = ['protocol_name']

    def __unicode__(self):
        return str(self.protocol_name)


class Seqtech(models.Model):
    seqtech_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=100)
    notes = models.TextField(blank=True)

class Lifestage(models.Model):
    lifestage_id = models.AutoField(primary_key=True)
    lifestage = models.CharField(unique=True, max_length=45)
    notes = models.CharField(max_length=400, default=None)
    class Meta:
        ordering = ['lifestage']
        verbose_name_plural="Lifecycle stages"

    def __unicode__(self):
        return str(self.lifestage)

class Phenotype(models.Model):
    phenotype_id = models.AutoField(primary_key=True)
    phenotype = models.CharField(unique=True, max_length=45)
    notes = models.CharField(max_length=45, default=None)
    def __unicode__(self):
        return unicode(self.phenotype)

class Growthphase(models.Model):
    growthphase = models.CharField(unique=True, max_length=100)
    notes = models.CharField(max_length=400, default=None, blank=True)

    class Meta:
        ordering = ['growthphase']

    def __unicode__(self):
        return unicode(self.growthphase)


class Genotype(models.Model):
    genotype = models.CharField(unique=True, max_length=45)
    notes = models.CharField(max_length=45, default=None, blank=True)

    class Meta:
        ordering = ['genotype']

    def __unicode__(self):
        return unicode(self.genotype)



class Collaborator(models.Model):
    collaborator_id = models.AutoField(primary_key=True)
    firstname = models.CharField(max_length=45)
    lastname = models.CharField(max_length=45)
    email = models.EmailField(unique=True, max_length=100)
    affiliation = models.CharField(max_length=100)
    sharepoint_site = models.URLField(blank=True)
    ftp_path = models.URLField(blank=True)
    ftp_username = models.CharField(max_length=50, blank=True)

    @property
    def name(self):
        return u"%s %s" % (self.firstname, self.lastname)

    def __unicode__(self):
        return unicode(self.name)


class Organism(models.Model):
    organism_id = models.AutoField(primary_key=True)
    organismcode = models.CharField(unique=True, max_length=10)
    genus = models.CharField(max_length=45)
    species = models.CharField(max_length=45)
    strain = models.CharField(max_length=45)
    def __unicode__(self):
        return unicode(self.organismcode)

class Software(models.Model):
    software_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    version = models.CharField(max_length=255)
    algorithm = models.CharField(max_length=255)
    source = models.CharField(max_length=255)
    sourceuri = models.CharField(max_length=255)
    notes = models.TextField(default=None)
    def __unicode__(self):
        return str(self.name)

class Analysistype(models.Model):
    analysistype_id = models.AutoField(primary_key=True)
    type = models.CharField(unique=True, max_length=255)
    definition = models.TextField()
    def __unicode__(self):
        return str(self.type)

class Dbxref(models.Model):
    dbxref_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=255)
    url = models.CharField(max_length=255)
    definition = models.TextField()
    def __unicode__(self):
        return str(self.dbxref_id)

class Cv(models.Model):
    cv_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=255, db_index=True)
    definition = models.TextField()
    def __unicode__(self):
        return str(self.name)

class Cvterm(models.Model):
    cvterm_id = models.AutoField(primary_key=True)
    cv = models.ForeignKey(Cv)
    name = models.CharField(max_length=1024, db_index=True)
    definition = models.TextField()
    dbxref = models.ForeignKey(Dbxref)
    is_obsolete = models.BooleanField(default="False")
    is_relationshiptype = models.BooleanField()
    def __unicode__(self):
        return str(self.name)

class Genome(models.Model):
    genome_id = models.AutoField(primary_key=True)
    organism = models.ForeignKey(Organism)
    version = models.CharField(max_length=45)
    source = models.CharField(max_length=45)
    svnrevision = models.CharField(max_length=45)
    svnpath = models.TextField(blank=True)
    def __unicode__(self):
        return str(self.organism)

class Result(models.Model):
    result_id = models.AutoField(primary_key=True)
    libraries = models.ManyToManyField('samples.Library')
    genome = models.ForeignKey(Genome)
    author = models.ForeignKey(Author)
    is_current = models.BooleanField(default="True")
    is_obsolete = models.BooleanField(default="False")
    analysispath = models.CharField(max_length=255)
    time_data_loaded = models.DateTimeField(auto_now=True)
    notes = models.TextField(default=None)
    def __unicode__(self):
        return str(self.result_id)

def get_resultfile_upload_destination(instance, filename):
    return "resultfiles/resultid_{id}/{file}".format(id=instance.result_id, file=filename)

class Result_CV(models.Model):
    result_cv_id = models.AutoField(primary_key=True)
    cvterm = models.TextField(db_index=True)
    definition = models.TextField()
    dbxref = models.ForeignKey(Dbxref)
    is_obsolete = models.BooleanField(default="False")
    is_relationshiptype = models.BooleanField()

class Resultfile(models.Model):
    resultfile_id = models.AutoField(primary_key=True)
    result = models.ForeignKey(Result)
    notes = models.CharField(max_length=1000, default='result')
    file = models.FileField(upload_to=get_resultfile_upload_destination)
    def __unicode__(self):
        return str(self.result)

class Resultprop(models.Model):
    resultprop_id = models.AutoField(primary_key=True)
    result = models.ForeignKey(Result)
    cvterm = models.ForeignKey(Cvterm)
    value = models.TextField()
    def __unicode__(self):
        return str(self.resultprop_id)

class Resultraw(models.Model):
    resultraw_id = models.AutoField(primary_key=True)
    result = models.ForeignKey(Result)
    chromosome = models.CharField(max_length=250, db_index=True)
    position = models.IntegerField(db_index=True)
    totalcount = models.IntegerField()
    posstrandcount = models.IntegerField()
    negstrandcount = models.IntegerField()
    majorstrand = models.IntegerField()
    time_data_loaded = models.DateTimeField(auto_now=True)
    def __unicode__(self):
        return str(self.resultraw_id)

class Resultslsite(models.Model):
    resultslsite_id = models.AutoField(primary_key=True)
    result = models.ForeignKey(Result)
    slid = models.CharField(max_length=45)
    chromosome = models.CharField(max_length=15, db_index=True)
    position = models.IntegerField(db_index=True)
    readcount = models.IntegerField()
    readstrand = models.CharField(max_length=5)
    geneid = models.CharField(max_length=45, db_index=True)
    genestrand = models.IntegerField()
    cdsstart = models.IntegerField()
    cdsend = models.IntegerField()
    putative5utr = models.IntegerField()
    rank = models.FloatField(db_index=True)
    intervallength = models.IntegerField()
    slpercent = models.DecimalField(max_digits=10, decimal_places=7)
    time_data_loaded = models.DateTimeField(auto_now=True)
    #geneid_current = models.CharField(max_length=45, db_index=True)

    def __unicode__(self):
        return str(self.resultslsite_id)

class Resultslgene(models.Model):
    resultslgene_id = models.AutoField(primary_key=True)
    result = models.ForeignKey(Result)
    geneid = models.CharField(max_length=45, db_index=True)
    genestrand = models.IntegerField()
    chromosome = models.CharField(max_length=25)
    cdsstart = models.IntegerField()
    cdsend = models.IntegerField()
    putative5utr = models.IntegerField()
    sensesitecount = models.IntegerField()
    antisensesitecount = models.IntegerField()
    sensereadcount = models.IntegerField()
    antisensereadcount = models.IntegerField()
    senseupsitecount = models.IntegerField()
    sensedownsitecount = models.IntegerField()
    time_data_loaded = models.DateTimeField(auto_now=True)
    #geneid_current = models.CharField(max_length=45, db_index=True)

    def __unicode__(self):
        return str(self.resultslgene_id)


class Resultsriboprof(models.Model):
    resultsriboprof_id = models.AutoField(primary_key=True)
    result = models.ForeignKey(Result)
    geneid = models.CharField(max_length=45, db_index=True)
    featuretype = models.CharField(max_length=10, db_index=True)
    counts_raw = models.DecimalField(max_digits=10, decimal_places=4)
    counts_normalized = models.DecimalField(max_digits=10, decimal_places=4)
    time_data_loaded = models.DateTimeField(auto_now=True)

    def __unicode__(self):
        return str(self.resultsriboprof_id)


class Analysis(models.Model):
    analysis_id = models.AutoField(primary_key=True)
    analysistype = models.ForeignKey(Analysistype)
    software = models.ForeignKey(Software)
    result = models.ForeignKey(Result)
    #libraries = models.ManyToManyField(Library) let's add this when we use South
    ordinal = models.IntegerField()
    time_data_loaded = models.DateTimeField(auto_now=True)
    notes = models.TextField(null=True)
    def __unicode__(self):
        return str(self.analysis_id)

def get_analysisfile_upload_destination(instance, filename):
    return "analysesfiles/resultid_{id}/{file}".format(id=instance.analysis.result_id, file=filename)


class AnalysisCV(models.Model):
	analysis_cv_id = models.AutoField(primary_key=True)
	cvterm = models.TextField(db_index=True)
	definition = models.TextField()
	dbxref = models.ForeignKey(Dbxref)
	is_obsolete = models.BooleanField(default="False")
	is_relationshiptype = models.BooleanField()

class Analysisfile(models.Model):
    analysisfile_id = models.AutoField(primary_key=True)
    analysis = models.ForeignKey(Analysis)
    notes = models.CharField(max_length=1000, default='analysis')
    file = models.FileField(upload_to=get_analysisfile_upload_destination)
    def __unicode__(self):
        return str(self.analysis)


class Analysisprop(models.Model):
    analysisprop_id = models.AutoField(primary_key=True)
    analysis = models.ForeignKey(Analysis)
    cvterm = models.ForeignKey(Cvterm)
    value = models.TextField()
    def __unicode__(self):
        return str(self.analysisprop_id)

class Feature(models.Model):
    feature_id = models.AutoField(primary_key=True)
    featuretype = models.CharField(max_length=40)
    genome = models.ForeignKey(Genome)
    chromosome = models.CharField(max_length=45, db_index=True)
    geneid = models.CharField(max_length=100, db_index=True)
    fmin = models.IntegerField()
    fmax = models.IntegerField()
    strand = models.IntegerField()
    phase = models.CharField(max_length=5)
    geneproduct = models.CharField(max_length=500, db_index=True)
    annotation = models.CharField(max_length=500, db_index=True)
    na_seq = models.TextField()
    aa_seq = models.TextField()
    time_data_loaded = models.DateTimeField(auto_now_add=True)
    time_data_modified = models.DateTimeField(auto_now=True)
    def __unicode__(self):
        return str(self.geneid)


class Geneidmap(models.Model):
    geneidmap_id = models.AutoField(primary_key=True)
    geneid_current = models.CharField(max_length=50, db_index=True)
    geneid_previous = models.CharField(max_length=50, db_index=True)
    db_soruce_previous = models.CharField(max_length=50)
    db_version_previous = models.CharField(max_length=50)
    db_soruce_current = models.CharField(max_length=50)
    db_version_current = models.CharField(max_length=50)
    time_data_loaded = models.DateTimeField(auto_now=True)


class Experiment(models.Model):
    experiment_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=25, db_index=True)
    type = models.CharField(max_length=25, choices=EXPERIMENT_TYPE_CHOICES)
    description = models.CharField(max_length=100)
    notes = models.TextField()
    libraries = models.ManyToManyField('samples.Library')
    date_created = models.DateTimeField(auto_now_add=True)
    date_modified = models.DateTimeField(auto_now=True)
    author_modified = models.ForeignKey(User)

    def __unicode__(self):
        return str(self.name)

#==============================================
#==============================================
# A temp table to change result_libraries from ngsdbview.library to samples.library
# This table and data are saved to refer back in future for any inconsistencies that may arise.

class Tempmtom(models.Model):
    result_id = models.IntegerField(db_index=True, unique=True)
    ngsdbview_libid = models.IntegerField()
    ngsdbview_libcode = models.CharField(max_length=10, db_index=True)
    samples_libid   = models.IntegerField()


class Templibraryprop(models.Model):
    library_id = models.IntegerField()
    cvterm_id = models.IntegerField()
    value = models.TextField()

#==============================================
#==============================================



from django.db import models
from django.contrib.auth.models import User

class UserProfile(models.Model):
    """This class extends User to add new fields and function definitions"""
    # This field is required.
    id = models.AutoField(primary_key=True)
    user = models.OneToOneField(User)
    libraries = models.ManyToManyField('samples.Library',default=['AH006'])
    def __unicode__(self):
        return str(self.user)

    #This is used to create a handle for a user to his or her profile
def create_user_profile(sender, instance, created, **kwargs):
    if created:
        UserProfile.objects.create(user=instance)



#========================================================
#========================================================
#========================================================


# after syncdb following manual changes were made to the database
#ALTER TABLE ngsdbview_geneidmap ALTER time_data_loaded SET DEFAULT now()
#ALTER TABLE ngsdbview_resultslgene ALTER time_data_loaded SET DEFAULT now()
#ALTER TABLE ngsdbview_resultslsite ALTER time_data_loaded SET DEFAULT now()
#ALTER TABLE ngsdbview_result ALTER time_data_loaded SET DEFAULT now()
#ALTER TABLE ngsdbview_analysis ALTER time_data_loaded SET DEFAULT now()
#ALTER TABLE ngsdbview_analysis ALTER na_seq SET NULL
#ALTER TABLE ngsdbview_analysis ALTER aa_seq SET NULL
#ALTER TABLE ngsdbview_analysis ALTER annotation SET NULL

#load initial data from '/Users/gramasamy/djcode/ngsdbview03/ngdbview_initialdata.sql'
#load genome data like this
#-- copy ngsdbview_feature (featuretype, genome_id, chromosome, geneid, fmin, fmax, strand, phase, geneproduct)from '/tmp/Linfantum_TriTrypDB-4.0_forNGSDB_featuretable_genomeid5.tab'
#-- copy ngsdbview_feature (featuretype, genome_id, chromosome, geneid, fmin, fmax, strand, phase, geneproduct)from '/tmp/Linfantum_TriTrypDB-2.0_forNGSDB_featuretable_genomeid1.tab'



