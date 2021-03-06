from snpdb.models import *
from samples.models import *
from django.shortcuts import render_to_response
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from utils import build_orderby_urls
from django.db.models import *
from django_boolean_sum import BooleanSum
from templatetags.snp_filters import *
from django.template import RequestContext
from GChartWrapper import *
from collections import *
import subprocess
import datetime
import os
import csv
import vcf
from operator import itemgetter


low_effects = ["SYNONYMOUS_START", "NON_SYNONYMOUS_START", "START_GAINED", "SYNONYMOUS_CODING", "SYNONYMOUS_STOP"]
high_effects = ["SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR", "START_LOST", "EXON_DELETED", "FRAME_SHIFT", "STOP_GAINED", "STOP_LOST", "RARE_AMINO_ACI"]
moderate_effects = ["NON_SYNONYMOUS_CODING", "CODON_CHANGE", "CODON_INSERTION", "CODON_CHANGE_PLUS_CODON_INSERTION",
                    "CODON_DELETION", "CODON_CHANGE_PLUS_CODON_DELETION", "UTR_5_DELETED", "UTR_3_DELETED"]
modifier_effects = ["UTR_5_PRIME", "UTR_3_PRIME", "REGULATION", "UPSTREAM", "DOWNSTREAM", "GENE", "TRANSCRIPT", "EXON",
                    "INTRON_CONSERVED", "INTRON", "INTRAGENIC", "INTERGENIC", "INTERGENIC_CONSERVED", "NONE", "CHROMOSOME", "CUSTOM", "CDS"]


def dashboard(request):
	title = "SNP Dashboard"

	lib_count = SNP.objects.values("library__library_code").distinct().annotate(Count('snp_id'))
	lib_snps = []
	lib_snp_total = 0
	for each in lib_count.iterator():
		lib_snps.append(each['snp_id__count'])
		lib_snp_total += each['snp_id__count']

	org_count = SNP.objects.values("result__genome__organism__organismcode").distinct().annotate(Count('snp_id'))
	org_snps = []
	org_snp_total = 0
	for each in org_count.iterator():
		org_snps.append(each['snp_id__count'])
		org_snp_total += each['snp_id__count']

	path = os.path.abspath(os.path.dirname(__file__))
	chart_path = os.path.join(path, 'gcharts/%s_impact.csv')
	image_path = 'snpdb/static/snps_by_%s.png'
	images_path = 'snps_by_%s.png'
	if os.path.isfile(chart_path % 'high') and os.path.isfile(image_path % 'high'):
		print "file ", chart_path % 'high', "and file ", image_path % 'high', " was found."
		pass
	else:
		save_snp_dashboard_files(chart_path, image_path)
	totals = [lib_snp_total, org_snp_total]

	#read count files
	high_count = read(chart_path % 'high')
	low_count = read(chart_path % 'low')
	moderate_count = read(chart_path % 'moderate')
	modifier_count = read(chart_path % 'modifier')
	impact_count = read(chart_path % 'impact')


	images = [images_path % 'library', images_path % 'organism', images_path % 'impact', images_path % 'high',
	          images_path % 'low', images_path % 'moderate', images_path % 'modifier']

	return render_to_response('snpdb/dashboard.html', {"title": title,
	                                                   "images": images,
	                                                   "totals": totals,
	                                                   "lib_count": lib_count,
	                                                   "org_count": org_count,
	                                                   "impact_count": impact_count,
	                                                   "high_count": high_count,
	                                                   "low_count": low_count,
	                                                   "moderate_count": moderate_count,
	                                                   "modifier_count": modifier_count,
	                                                   },  context_instance=RequestContext(request))


# Returns the general effect table view.
def effect(request):
	order_by = request.GET.get('order_by', 'effect')
	current_url = request.get_full_path()
	snp_effect_list = Effect.objects.all().order_by(order_by)
	paginator = Paginator(snp_effect_list, 50)
	page = request.GET.get('page')

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['snp_id', 'effect', 'effect_class',
	                                                           'effect_string', 'effect_group'])
	try:
		snp_effect = paginator.page(page)
	except PageNotAnInteger:
		snp_effect = paginator.page(1)
	except EmptyPage:
		snp_effect = paginator.page(paginator.num_pages)

	toolbar_max = min(snp_effect.number + 4, paginator.num_pages)
	toolbar_min = max(snp_effect.number - 4, 0)

	return render_to_response('snpdb/effect.html', {"snp_effect": snp_effect,
	                                                "filter_urls": filter_urls,
	                                                "paginator": paginator,
	                                                "toolbar_max": toolbar_max,
	                                                "toolbar_min": toolbar_min,
	                                                "current_url": current_url},
	                          context_instance=RequestContext(request))


# Returns the general filter table view
def snp_filter(request):
	order_by = request.GET.get('order_by', 'snp')
	filter_list = Filter.objects.all().order_by(order_by)
	paginator = Paginator(filter_list, 50)

	page = request.GET.get('page')

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['snp_id', 'filter_id', 'filter_result',
	                                                           'filter_cv'])
	try:
		filters = paginator.page(page)
	except PageNotAnInteger:
		filters = paginator.page(1)
	except EmptyPage:
		filters = paginator.page(paginator.num_pages)

	toolbar_max = min(filters.number + 3, paginator.num_pages)
	toolbar_min = max(filters.number - 3, 0)
	return render_to_response('snpdb/filter.html', {"filters": filters,
	                                                "filter_urls": filter_urls,
	                                                "paginator": paginator,
	                                                "toolbar_max": toolbar_max,
	                                                "toolbar_min": toolbar_min})


# Returns the general SNP table view
def snp(request):
	order_by = request.GET.get('order_by', 'snp_id')
	snp_list = SNP.objects.values('snp_id', 'snp_position', 'result', 'ref_base', 'alt_base',
	                              'heterozygosity', 'quality', 'library__library_code', 'chromosome__chromosome_name').order_by(order_by)
	count = len(snp_list)
	paginator = Paginator(snp_list, 50)
	page = request.GET.get('page')

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['snp_id', 'snp_position', 'result',
	                                                           'ref_base', 'alt_base', 'heterozygosity',
	                                                           'quality', 'library__library_code', 'chromosome__chromosome_name'])
	try:
		snps = paginator.page(page)
	except PageNotAnInteger:
		snps = paginator.page(1)
	except EmptyPage:
		snps = paginator.page(paginator.num_pages)

	toolbar_max = min(snps.number + 3, paginator.num_pages)
	toolbar_min = max(snps.number - 3, 0)

	return render_to_response('snpdb/snp.html', {"snps": snps,
	                                             "count": count,
	                                             "filter_urls": filter_urls,
	                                             "paginator": paginator,
	                                             "toolbar_max": toolbar_max,
	                                             "toolbar_min": toolbar_min})


#Lists CNV values
def cnv(request):
	order_by = request.GET.get('order_by', 'cnv_id')
	cnv_list = CNV.objects.all().order_by(order_by)
	count = len(cnv_list)
	paginator = Paginator(cnv_list, 50)
	page = request.GET.get('page')

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['cnv_id', 'chromosome__chromosome_name', 'coordinte',
	                                                           'CNV_value', 'result_id', 'library__library_code'])
	try:
		cnvs = paginator.page(page)
	except PageNotAnInteger:
		cnvs = paginator.page(1)
	except EmptyPage:
		cnvs = paginator.page(paginator.num_pages)

	toolbar_max = min(cnvs.number + 3, paginator.num_pages)
	toolbar_min = max(cnvs.number - 3, 0)

	return render_to_response('snpdb/CNV.html', {"cnvs": cnvs,
	                                             "count": count,
	                                             "filter_urls": filter_urls,
	                                             "paginator": paginator,
	                                             "toolbar_max": toolbar_max,
	                                             "toolbar_min": toolbar_min})

# Returns the general SNP Type table view.
def snp_type(request):
	order_by = request.GET.get('order_by', 'snptype_id')
	snptype_list = SNP_Type.objects.all().order_by(order_by)
	paginator = Paginator(snptype_list, 50)

	page = request.GET.get('page')

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['snptype_id', 'snp_id', 'indel',
	                                                           'deletion', 'is_snp', 'monomorphic',
	                                                           'transition', 'sv'])
	try:
		snptypes = paginator.page(page)
	except PageNotAnInteger:
		snptypes = paginator.page(1)
	except EmptyPage:
		snptypes = paginator.page(paginator.num_pages)

	toolbar_max = min(snptypes.number + 4, paginator.num_pages)
	toolbar_min = max(snptypes.number - 4, 0)

	return render_to_response('snpdb/snptype.html', {"snptypes": snptypes,
	                                                 "filter_urls": filter_urls,
	                                                 "paginator": paginator,
	                                                 "toolbar_max": toolbar_max,
	                                                 "toolbar_min": toolbar_min})


# Returns the general statistics table view.
def statistics(request):
	order_by = request.GET.get('order_by', 'stats_id')
	statistic_list = Statistics.objects.all().order_by(order_by)
	paginator = Paginator(statistic_list, 50)

	page = request.GET.get('page')

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['stats_id', 'snp',
	                                                           'stats_cvterm', 'cv_value'])
	try:
		statistic = paginator.page(page)
	except PageNotAnInteger:
		statistic = paginator.page(1)
	except EmptyPage:
		statistic = paginator.page(paginator.num_pages)

	toolbar_max = min(statistic.number + 3, paginator.num_pages)
	toolbar_min = max(statistic.number - 3, 0)

	return render_to_response('snpdb/statistics.html', {"statistic": statistic,
	                                                    "filter_urls": filter_urls,
	                                                    "paginator": paginator,
	                                                    "toolbar_max": toolbar_max,
	                                                    "toolbar_min": toolbar_min})


# Search views
#---------------------------------------------------------------------------------------------------
#todo merge filter views with basic views.
# Returns a view of the Effect table that has been filtered.
def effect_filter(request):
	selection = request.GET.get('att')
	filter_on = request.GET.get('s')
	# current_url = request.get_full_path()
	filter_dict = {}
	filter_dict[str(selection)] = str(filter_on)
	result_list = Effect.objects.all().filter(**filter_dict)

	order_by = request.GET.get('order_by', 'snp')
	result_list = result_list.order_by(order_by)
	paginator = Paginator(result_list, 50)
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['snp_id', 'effect', 'effect_class',
	                                                           'effect_string', 'effect_group'])
	try:
		snp_effect = paginator.page(page)
	except PageNotAnInteger:
		snp_effect = paginator.page(1)
	except EmptyPage:
		snp_effect = paginator.page(paginator.num_pages)

	toolbar_max = min(snp_effect.number + 3, paginator.num_pages)
	toolbar_min = max(snp_effect.number - 3, 0)

	return render_to_response('snpdb/effect.html', {"snp_effect": snp_effect,
	                                                "filter_urls": filter_urls,
	                                                "selection": selection,
	                                                "filter_on": filter_on,
	                                                "paginator": paginator,
	                                                "toolbar_max": toolbar_max,
	                                                "toolbar_min": toolbar_min},
	                          context_instance=RequestContext(request))


# Returns a view of the SNP table that has been filtered.
def snp_filter_result(request):
	selection = request.GET.get('att')
	filter_on = request.GET.get('s')
	filter_dict = {}
	filter_dict[str(selection)] = str(filter_on)

	result_list = SNP.objects.values('snp_id', 'snp_position', 'result', 'ref_base', 'alt_base',
	                                 'heterozygosity', 'quality', 'library__library_code', 'chromosome__chromosome_name').filter(**filter_dict)

	order_by = request.GET.get('order_by', 'snp_id')
	result_list = result_list.order_by(order_by)
	filter_urls = build_orderby_urls(request.get_full_path(), ['snp_id', 'snp_position', 'result',
	                                                           'ref_base', 'alt_base', 'heterozygosity',
	                                                           'quality', 'library__library_code', 'chromosome__chromosome_name'])
	paginator = Paginator(result_list, 50)
	page = request.GET.get('page')

	try:
		snps = paginator.page(page)
	except PageNotAnInteger:
		snps = paginator.page(1)
	except EmptyPage:
		snps = paginator.page(paginator.num_pages)

	toolbar_max = min(snps.number + 3, paginator.num_pages)
	toolbar_min = max(snps.number - 3, 0)

	return render_to_response('snpdb/snp.html', {"snps": snps,
	                                             "filter_urls": filter_urls,
	                                             "selection": selection,
	                                             "filter_on": filter_on,
	                                             "paginator": paginator,
	                                             "toolbar_max": toolbar_max,
	                                             "toolbar_min": toolbar_min},
	                          context_instance=RequestContext(request))


# Returns a view of the Filter table that has been filtered.
def filter_filter(request):
	selection = request.GET.get('att')
	filter_on = request.GET.get('s')
	filter_dict = {}
	filter_dict[str(selection)] = str(filter_on)
	result_list = Filter.objects.all().filter(**filter_dict)

	order_by = request.GET.get('order_by', 'snp')
	result_list = result_list.order_by(order_by)
	paginator = Paginator(result_list, 50)
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['snp_id', 'filter_id', 'filter_result',
	                                                           'filter_cv'])
	try:
		filters = paginator.page(page)
	except PageNotAnInteger:
		filters = paginator.page(1)
	except EmptyPage:
		filters = paginator.page(paginator.num_pages)

	toolbar_max = min(filters.number + 3, paginator.num_pages)
	toolbar_min = max(filters.number - 3, 0)

	return render_to_response('snpdb/filter.html', {"filters": filters,
	                                                "filter_urls": filter_urls,
	                                                "selection": selection,
	                                                "filter_on": filter_on,
	                                                "paginator": paginator,
	                                                "toolbar_max": toolbar_max,
	                                                "toolbar_min": toolbar_min})


# Returns a view of the SNP Type table that has been filtered.
def snptype_filter(request):
	selection = request.GET.get('att')
	filter_on = request.GET.get('s')
	filter_dict = {}
	filter_dict[str(selection)] = str(filter_on)

	result_list = SNP_Type.objects.all().filter(**filter_dict)

	order_by = request.GET.get('order_by', 'snp')
	result_list = result_list.order_by(order_by)
	paginator = Paginator(result_list, 50)
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['snptype_id', 'snp_id', 'indel',
	                                                           'deletion', 'is_snp', 'monomorphic',
	                                                           'transition', 'sv'])
	try:
		snptypes = paginator.page(page)
	except PageNotAnInteger:
		snptypes = paginator.page(1)
	except EmptyPage:
		snptypes = paginator.page(paginator.num_pages)

	toolbar_max = min(snptypes.number + 4, paginator.num_pages)
	toolbar_min = max(snptypes.number - 4, 0)

	return render_to_response('snpdb/snptype.html', {"snptypes": snptypes,
	                                                 "filter_urls": filter_urls,
	                                                 "selection": selection,
	                                                 "filter_on": filter_on,
	                                                 "paginator": paginator,
	                                                 "toolbar_max": toolbar_max,
	                                                 "toolbar_min": toolbar_min})


# Returns a view of the Statistics table that has been filtered.
def statistics_filter(request):
	selection = request.GET.get('att')
	filter_on = request.GET.get('s')
	filter_dict = {}
	filter_dict[str(selection)] = str(filter_on)

	result_list = Statistics.objects.all().filter(**filter_dict)

	order_by = request.GET.get('order_by', 'snp')
	result_list = result_list.order_by(order_by)
	paginator = Paginator(result_list, 50)
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['stats_id', 'snp', 'stats_cvterm',
	                                                           'cv_value'])
	try:
		statistic = paginator.page(page)
	except PageNotAnInteger:
		statistic = paginator.page(1)
	except EmptyPage:
		statistic = paginator.page(paginator.num_pages)

	toolbar_max = min(statistic.number + 4, paginator.num_pages)
	toolbar_min = max(statistic.number - 4, 0)

	return render_to_response('snpdb/statistics.html', {"statistic": statistic,
	                                                    "filter_urls": filter_urls,
	                                                    "selection": selection,
	                                                    "filter_on": filter_on,
	                                                    "paginator": paginator,
	                                                    "toolbar_max": toolbar_max,
	                                                    "toolbar_min": toolbar_min})


# Query views.
# --------------------------------------------------------------------------------------
# Displays the search page to compare snps across libraries
def compare_gene_lib(request):
	return render_to_response('snpdb/compare_gene_library.html', )


# Returns a list of libraries that the desired gene is found in.
def compare_gene_lib_filter(request):
	order_by = request.GET.get('order_by', 'library__library_code')
	gene = request.GET.get('s')
	genes = gene.split()

	genomes = Feature.objects.values_list('genome_id', flat=True).filter(geneid__in=genes).distinct()
	libraries = []
	for each in genomes:
		result_list = SNP.objects.values('library__library_code', 'result__genome__organism__organismcode', 'result__genome__version', 'result__genome__genome_id').filter(result__genome__genome_id=each).distinct().order_by(order_by)
		libraries.append(result_list)
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['library__library_code',
	                                                           'result__genome__organism__organismcode',
	                                                           'result__genome____genome_id', 'result__genome__version'])
	count = len(libraries[0])
	paginator = Paginator(libraries, 50)
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/compare_gene_library_filter.html', {"results": results,
	                                                                     "gene": gene,
	                                                                     "count": count,
	                                                                     "filter_urls": filter_urls,
	                                                                     "toolbar_max": toolbar_max,
	                                                                     "toolbar_min": toolbar_min})


# Returns the comparison of a gene across specific libraries.
def compare_gene_lib_filter_results_effect(request):
	# order_by = request.GET.get('order_by', 'library__library_code')
	gene_string = request.GET.get('s')
	genes = gene_string.split()
	library = request.GET.getlist('check')
	test = {}
	for gene in genes:
		cds_fmin = Feature.objects.values_list('fmin', flat=True).filter(geneid=gene, featuretype='CDS')[0]
		cds_fmax = Feature.objects.values_list('fmax', flat=True).filter(geneid=gene, featuretype='CDS')[0]
		fmin = Feature.objects.filter(geneid=gene).filter(featuretype='gene').values('fmin')[0]
		fmax = Feature.objects.filter(geneid=gene).filter(featuretype='gene').values('fmax')[0]
		chromosome = Feature.objects.filter(geneid=gene).filter(featuretype='gene').values_list('chromosome', flat=True)[0]

		result_list = SNP.objects.filter(effect__effect_id=6, effect__effect_string__exact=gene,
		                                 library__library_code__in=library).values('library', 'library__library_code', 'snp_id',
		                                                                           'snp_position', 'ref_base', 'alt_base',
		                                                                           'heterozygosity', 'quality',
		                                                                           'chromosome__chromosome_name', 'effect__effect_string',
		                                                                           'effect__effect_class', 'effect__effect').distinct().order_by('snp_position')
		#Checks to see if tuples have all libraries present. Inserts blank tuples if not.
		for each in result_list:
			new_tuple = [(None, None, None, gene, cds_fmin, cds_fmax, fmin['fmin'], fmax['fmax'], chromosome, each['snp_position'])] * len(library)
			curr_library = each['library__library_code']
			tup = (curr_library, each['ref_base'], each['alt_base'], gene, cds_fmin, cds_fmax, fmin['fmin'], fmax['fmax'], chromosome, each['snp_position'])
			index = library.index(curr_library)
			if each['snp_position'] in test:
				current_tup = test[each['snp_position']]
				current_tup[index] = tup
				test[each['snp_position']] = current_tup
			else:
				new_tuple[index] = tup
				test[each['snp_position']] = new_tuple
	test = OrderedDict(sorted(test.items(), key=lambda key: key[0]))
	count = len(test)
	filter_urls = build_orderby_urls(request.get_full_path(), ['gene', 'snp_position', 'ref_base', 'alt_base', 'library__library_code', 'fmin', 'fmax'])

	paginator = Paginator(test.items(), 50)
	page = request.GET.get('page')
	try:
		results = paginator.page(page)
		print "got this page"
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)
	return render_to_response('snpdb/compare_gene_library_filter_result.html', {"results": results,
	                                                                            "genes": genes,
	                                                                            "gene": gene_string,
	                                                                            "library": library,
	                                                                            "test": test,
	                                                                            "count": count,
	                                                                            "library_group": sorted(library),
	                                                                            "filter_urls": filter_urls,
	                                                                            "toolbar_max": toolbar_max,
	                                                                            "toolbar_min": toolbar_min}, context_instance=RequestContext(request))


# Returns the list of genes found within the selected libraries.
# todo need to change so that it references effect table rather than feature table
def compare_gene_lib_filter_results(request):
	order_by = request.GET.get('order_by', 'library__library_code')
	gene = request.GET.get('s')
	print gene
	library = request.GET.getlist('check')

	#Gets the start and stop position of the coding region of the gene.
	cds_fmin = Feature.objects.values_list('fmin', flat=True).filter(geneid=gene, featuretype='CDS')[0]
	cds_fmax = Feature.objects.values_list('fmax', flat=True).filter(geneid=gene, featuretype='CDS')[0]

	#Gets the start and stop position of gene.
	fmin = Feature.objects.filter(geneid=gene).filter(featuretype='gene').values('fmin')[0]
	fmax = Feature.objects.filter(geneid=gene).filter(featuretype='gene').values('fmax')[0]

	#Collects libraries that are effected by this gene
	result_list = SNP.objects.values('snp_id', 'snp_position',
	                                 'ref_base', 'alt_base',
	                                 'library__library_code').filter(library__library_code__in=library,
	                                                                 snp_position__range=(Feature.objects.values_list('fmin', flat=True).filter(geneid=gene).filter(featuretype='gene')[0],
	                                                                                      Feature.objects.values_list('fmax', flat=True).filter(geneid=gene).filter(featuretype='gene')[0]),
	                                                                 chromosome__chromosome_name=Feature.objects.values_list('chromosome', flat=True).filter(geneid=gene).filter(featuretype='gene')[0]).order_by(order_by)
	snp_group = []
	library_group = []
	max_snps = 0
	count = 0
	for each in result_list:
		snp_position = each['snp_position']
		library = each['library__library_code']
		if snp_position in snp_group:
			count += 1
			pass
		else:
			snp_group.append(snp_position)
			count = 1
			pass
		if count > max_snps:
			max_snps = count
			count = 0
		if library in library_group:
			pass
		else:
			library_group.append(library)
	paginator = Paginator(result_list, 50)
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['snp_id', 'snp_position', 'ref_base', 'alt_base', 'library__library_code'])
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/compare_gene_library_filter_result.html', {"results": results,
	                                                                            "gene": gene,
	                                                                            "cds_fmin": cds_fmin,
	                                                                            "cds_fmax": cds_fmax,
	                                                                            "fmin": fmin,
	                                                                            "fmax": fmax,
	                                                                            "count": len(snp_group),
	                                                                            "library_group": library_group,
	                                                                            "snp_group": snp_group,
	                                                                            "max_snps": range(max_snps),
	                                                                            "filter_urls": filter_urls,
	                                                                            "toolbar_max": toolbar_max,
	                                                                            "toolbar_min": toolbar_min})


# Returns information about a gene through the feature table.
def gene_feature(request):
	geneid = request.GET.get('geneid')
	order_by = request.GET.get('order_by', 'geneid')

	feature = Feature.objects.all().filter(geneid=geneid, featuretype='gene').order_by(order_by)
	paginator = Paginator(feature, 50)
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['snp_id', 'snp_position', 'ref_base',
	                                                           'alt_base', 'library__library_code'])
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/gene_feature.html', {"geneid": geneid,
	                                                      "results": results,
	                                                      "filter_urls": filter_urls,
	                                                      "toolbar_max": toolbar_max,
	                                                      "toolbar_min": toolbar_min,})


# The search view for the user to input a gene. Lists all gene ids for the user to choose from.
def gene_snps(request):
	# genes = Effect.objects.values('effect_string').filter(effect=6).filter(effect_class=("NON_SYNONYMOUS_CODING" or "SYNONYMOUS_CODING")).distinct().order_by('effect_string')
	# paginator = Paginator(genes, 120)

	# page = request.GET.get('page')
	# try:
	# 	genes = paginator.page(page)
	# except PageNotAnInteger:
	# 	genes = paginator.page(1)
	# except EmptyPage:
	# 	genes = paginator.page(paginator.num_pages)
	# toolbar_max = min(genes.number + 4, paginator.num_pages)
	# toolbar_min = max(genes.number - 4, 0)

	return render_to_response('snpdb/gene_to_snp.html')


# Returns all snps found within the gene location regardless of library.
def gene_snps_filter(request):
	flanks = int(request.GET.get('f'))
	order_by = request.GET.get('order_by', 'library__library_code')
	gene = request.GET.get('s')
	fmin = Feature.objects.values_list('fmin', flat=True).filter(geneid=gene).filter(featuretype='gene')[0]
	fmax = Feature.objects.values_list('fmax', flat=True).filter(geneid=gene).filter(featuretype='gene')[0]
	try:
		cds_fmin = Feature.objects.values_list('fmin', flat=True).filter(geneid=gene, featuretype='CDS')[0]
		cds_fmax = Feature.objects.values_list('fmax', flat=True).filter(geneid=gene, featuretype='CDS')[0]
	except IndexError:
		print "There are no CDS regions associated with this gene."
		cds_fmin = 0
		cds_fmax = 0
		pass
	# chromosome = Feature.objects.values_list('chromosome', flat=True).filter(geneid=gene).filter(featuretype='gene')[0]
	result_list = SNP.objects.values('library__library_code', 'result_id',
	                                 'chromosome__chromosome_name', 'snp_id',
	                                 'snp_position', 'ref_base',
	                                 'alt_base').filter(snp_position__range=((Feature.objects.values_list('fmin', flat=True).filter(geneid=gene).filter(featuretype='gene')[0])+flanks,
	                                                                         (Feature.objects.values_list('fmax', flat=True).filter(geneid=gene).filter(featuretype='gene')[0])+flanks),
	                                                    chromosome__chromosome_name=Feature.objects.values_list('chromosome',
	                                                                                                            flat=True).filter(geneid=gene).filter(featuretype='gene')[0]).order_by(order_by)
	count = result_list.count()
	paginator = Paginator(result_list, 50)
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['library__library_code', 'result_id',
	                                                           'chromosome__chromosome_name', 'snp_id',
	                                                           'snp_position', 'ref_base', 'alt_base'])
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/snps_in_gene_filter.html', {"results": results,
	                                                             "filter_urls": filter_urls,
	                                                             "paginator": paginator,
	                                                             "toolbar_max": toolbar_max,
	                                                             "toolbar_min": toolbar_min,
	                                                             "cds_fmin": cds_fmin,
	                                                             "cds_fmax": cds_fmax,
	                                                             "fmin": fmin,
	                                                             "fmax": fmax,
	                                                             "gene": gene,
	                                                             "count": count})


# Returns the search page to query both a library/gene combination for snps.
def library_gene_snps(request):
	gene = request.GET.get('s')
	library = request.GET.get('lib')
	lib_list = Library.objects.values('library_code').order_by('library_code')
	page = request.GET.get('page')
	paginator = Paginator(lib_list, 120)

	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)
	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)
	return render_to_response('snpdb/library_to_snp.html', {"results": results,
	                                                        "gene": gene,
	                                                        "library": library,
	                                                        "paginator": paginator,
	                                                        "toolbar_max": toolbar_max,
	                                                        "toolbar_min": toolbar_min})


# Returns the snps found within a specific library and gene.
def library_gene_snps_filter(request):
	order_by = request.GET.get('order_by', 'library__library_code')
	gene = request.GET.get('s')
	library = request.GET.get('lib')
	cds_fmin = Feature.objects.values_list('fmin', flat=True).filter(geneid=gene, featuretype='CDS')[0]
	cds_fmax = Feature.objects.values_list('fmax', flat=True).filter(geneid=gene, featuretype='CDS')[0]
	fmin = Feature.objects.values_list('fmin', flat=True).filter(geneid=gene).filter(featuretype='gene')[0]
	fmax = Feature.objects.values_list('fmax', flat=True).filter(geneid=gene).filter(featuretype='gene')[0]
	result_list = SNP.objects.values('library__library_code', 'result_id',
	                                 'chromosome__chromosome_name', 'snp_id',
	                                 'snp_position', 'ref_base',
	                                 'alt_base', 'quality',
	                                 'heterozygosity').filter(snp_position__range=(Feature.objects.values_list('fmin', flat=True).filter(geneid=gene).filter(featuretype='gene')[0],
	                                                                               Feature.objects.values_list('fmax', flat=True).filter(geneid=gene).filter(featuretype='gene')[0]),
	                                                          library__library_code=library,
	                                                          chromosome__chromosome_name=Feature.objects.values_list('chromosome',
	                                                                                                                  flat=True).filter(geneid=gene).filter(featuretype='gene')[0]).order_by(order_by)
	count = result_list.count()
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['library__library_code', 'result_id',
	                                                           'chromosome__chromosome_name', 'snp_id',
	                                                           'snp_position', 'ref_base', 'alt_base',
	                                                           'quality', 'heterozygosity'])
	paginator = Paginator(result_list, 50)
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)
		print "Error"

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	print toolbar_max
	return render_to_response('snpdb/library_to_snp_filter.html', {"results": results,
	                                                               "gene": gene,
	                                                               "library": library,
	                                                               "cds_fmin": cds_fmin,
	                                                               "cds_fmax": cds_fmax,
	                                                               "fmin": fmin,
	                                                               "fmax": fmax,
	                                                               "filter_urls": filter_urls,
	                                                               "paginator": paginator,
	                                                               "toolbar_max": toolbar_max,
	                                                               "toolbar_min": toolbar_min,
	                                                               "count": count})


# Returns a summary of the number of snps found in each library.
def library_snp_summary(request):
	order_by = request.GET.get('order_by', 'library')
	results = SNP.objects.values('library_id',
	                             'library__library_code', 'result__genome__organism__organismcode',
	                             ).distinct().annotate(num_snps=Count('snp_id'),
	                                                   hetero=BooleanSum('heterozygosity'),
	                                                   indel=BooleanSum('snp_type__indel'),
	                                                   trans=BooleanSum('snp_type__transition'),
	                                                   ).order_by(order_by)
	result = []
	for each in results:
		library=each['library__library_code']
		organismcode = SNP.objects.values_list('result__genome__organism__organismcode', flat=True).filter(library__library_code=library)[0]
		contig = get_chromosome_size(organismcode)
		for key, value in contig.iteritems():
			each[key] = value
		result.append(each)
	paginator = Paginator(result, 50)
	page = request.GET.get('page')
	count = len(result)

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['library__library_id', 'library__library_code', 'snp_id', 'num_snps', 'hetero', 'homo', 'indel', 'trans', 'snp_density'])
	try:
		result = paginator.page(page)
	except PageNotAnInteger:
		result = paginator.page(1)
	except EmptyPage:
		result = paginator.page(paginator.num_pages)

	toolbar_max = min(result.number + 4, paginator.num_pages)
	toolbar_min = max(result.number - 4, 0)

	return render_to_response('snpdb/library_snp_summary.html', {"result": result,
	                                                             "count": count,
	                                                             "filter_urls": filter_urls,
	                                                             "paginator": paginator,
	                                                             "toolbar_max": toolbar_max,
	                                                             "toolbar_min": toolbar_min})


# Returns all snps found within a library using the effect table.
def library_snps(request):
	order_by = request.GET.get('order_by', 'snp_id').encode("ascii")
	library = request.GET.get('lib')
	count = request.GET.get('count')
	selection = request.GET.get('att')
	filter_on = request.GET.get('s')
	filter_dict = {}
	if selection:
		filter_dict[str(selection)] = str(filter_on)
	if filter_dict:
		if selection == 'effect__effect_string':
			results = SNP.objects.filter(effect__effect_id=6, effect__effect_string__exact=filter_on.decode('utf-8'),
			                             effect__effect_class__endswith='SYNONYMOUS_CODING'.decode('utf-8'),
			                             library__library_code=library).values('library', 'library__library_code', 'snp_id',
			                                                                   'snp_position', 'ref_base', 'alt_base',
			                                                                   'heterozygosity', 'quality',
			                                                                   'chromosome__chromosome_name', 'effect__effect_string',
			                                                                   'effect__effect_class', 'effect__effect')
		else:
			results = SNP.objects.filter(**filter_dict).filter(library__library_code=library).values('library', 'library__library_code', 'snp_id',
			                                                                                         'snp_position', 'ref_base', 'alt_base',
			                                                                                         'heterozygosity', 'quality',
			                                                                                         'chromosome__chromosome_name', 'effect__effect_string',
			                                                                                         'effect__effect_class', 'effect__effect')
	else:
		results = SNP.objects.values('library', 'library__library_code', 'snp_id',
		                             'snp_position', 'ref_base', 'alt_base',
		                             'heterozygosity', 'quality',
		                             'chromosome__chromosome_name', 'effect__effect_string',
		                             'effect__effect_class', 'effect__effect')

	# snp_dict = {}
	# for each in results:
	#     current_genes = []
	#     if each['library__library_code'] == library:
	#         if empty_effect(each['effect__effect']) or each['effect__effect'] == 6:
	#             if empty_effect(each['effect__effect_class']) or (each['effect__effect_class'] == ('NON_SYNONYMOUS_CODING' or 'SYNONYMOUS_CODING')):
	#                 if each['snp_id'] in snp_dict:
	#                     for k, v in snp_dict[each['snp_id']].iteritems():
	#                         if k == 'effect__effect_string':
	#                             if type(v) is list:
	#                                 for x in v:
	#                                     current_genes.append(str(x).decode('UTF8').strip())
	#                             else:
	#                                 current_genes.append(str(v).decode('UTF8'))
	#                     if each['effect__effect_string'] in current_genes:
	#                         pass
	#                     elif snp_dict[each['snp_id']]['effect__effect_string'] == 'None':
	#                         snp_dict[each['snp_id']] = each
	#                     else:
	#                         current_genes.append(str(each['effect__effect_string'].decode('UTF8').strip()))
	#                         each['effect__effect_string'] = current_genes
	#                         snp_dict[each['snp_id']] = each
	#                 else:
	#                     snp_dict[each['snp_id']] = each
	#             else:
	#                 if each['snp_id'] in snp_dict:
	#                     pass
	#                 else:
	#                     each["effect__effect_class"] = 'None'
	#                     each["effect__effect_string"] = 'None'
	#                     each["effect__effect"] = 'None'
	#                     snp_dict[each['snp_id']] = each
	#         else:
	#             if each['snp_id'] in snp_dict:
	#                 pass
	#             else:
	#                 each["effect__effect_class"] = 'None'
	#                 each["effect__effect_string"] = 'None'
	#                 each["effect__effect"] = 'None'
	#                 snp_dict[each['snp_id']] = each
	# sorted_snp_dict = sorted(snp_dict.items(), key=lambda y: y[1][order_by])
	sorted_snp_dict = genes_from_effect(results, library, order_by)
	paginator = Paginator(sorted_snp_dict, 100)
	page = request.GET.get('page')

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['library', 'library__library_code', 'snp_id',
	                                                           'snp_position', 'ref_base', 'alt_base',
	                                                           'heterozygosity', 'quality',
	                                                           'chromosome__chromosome_name',
	                                                           'effect__effect_string'])
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/library_snps.html', {"results": results,
	                                                      "library": library,
	                                                      "order_by": order_by,
	                                                      "filter_urls": filter_urls,
	                                                      "paginator": paginator,
	                                                      "toolbar_max": toolbar_max,
	                                                      "toolbar_min": toolbar_min,
	                                                      "count": count})


def genes_from_effect(results, library, order_by):
	snp_dict = {}
	for each in results:
		# print each
		current_genes = []
		if each['library__library_code'] == library:
			if empty_effect(each['effect__effect']) or each['effect__effect'] == 6:
				# if empty_effect(each['effect__effect_class']) or (each['effect__effect_class'] == ('NON_SYNONYMOUS_CODING' or 'SYNONYMOUS_CODING')):
				if each['snp_id'] in snp_dict:
					for k, v in snp_dict[each['snp_id']].iteritems():
						if k == 'effect__effect_string':
							if type(v) is list:
								for x in v:
									current_genes.append(str(x).decode('UTF8').strip())
							else:
								current_genes.append(str(v).decode('UTF8'))
					if each['effect__effect_string'] in current_genes:
						pass
					elif snp_dict[each['snp_id']]['effect__effect_string'] == 'None':
						snp_dict[each['snp_id']] = each
					else:
						current_genes.append(str(each['effect__effect_string']).strip())
						each['effect__effect_string'] = current_genes
						snp_dict[each['snp_id']] = each
				else:
					snp_dict[each['snp_id']] = each
				# else:
				#     if each['snp_id'] in snp_dict:
				#         pass
				#     else:
				#         each["effect__effect_class"] = 'None'
				#         each["effect__effect_string"] = 'None'
				#         each["effect__effect"] = 'None'
				#         snp_dict[each['snp_id']] = each
			else:
				if each['snp_id'] in snp_dict:
					pass
				else:
					each["effect__effect_class"] = 'None'
					each["effect__effect_string"] = 'None'
					each["effect__effect"] = 'None'
					snp_dict[each['snp_id']] = each
	sorted_snp_dict = sorted(snp_dict.items(), key=lambda x: x[1][order_by])
	return sorted_snp_dict


#todo determine if snps are in coding region of gene.
# Returns snps found in a library and chromosome.
def library_chromosome_snps_filter(request):
	chromosome = request.GET.get('s')
	library = request.GET.get('lib')
	order_by = request.GET.get('order_by', 'snp_position')
	genome = SNP.objects.values_list('result__genome__genome_id', flat=True).filter(library__library_code=library).distinct()[0]

	#Collects the ranges of all genes for the specific chromosome
	ranges = Feature.objects.values_list('fmin', 'fmax', 'chromosome').filter(chromosome=chromosome).filter(featuretype='gene', genome_id=genome).order_by('fmin')
	# result_list = []
	for each in ranges:
		#Collects snps found within the gene range
		# results = SNP.objects.values('library__library_code', 'result_id',
		#                              'chromosome__chromosome_name',
		#                              'snp_position', 'ref_base', 'alt_base', 'quality',
		#                              'heterozygosity').filter(snp_position__range=(each[0],
		#                                                                            each[1]),
		#                                                       library__library_code=library,
		#                                                       chromosome__chromosome_name=each[2]).order_by(order_by)
		# if results is not None:
		# 	result_list.append(results)

		#Collects all genes found within the chromosome
		results = SNP.objects.values('library__library_code', 'result_id',
		                             'chromosome__chromosome_name',
		                             'snp_position', 'ref_base', 'alt_base', 'quality',
		                             'heterozygosity').filter(library__library_code=library,
		                                                      chromosome__chromosome_name=each[2]).order_by(order_by)

	count = len(results)
	page = request.GET.get('page')
	filter_urls = build_orderby_urls(request.get_full_path(), ['library__library_code', 'result_id',
	                                                           'chromosome__chromosome_name', 'snp_id',
	                                                           'snp_position', 'ref_base', 'alt_base',
	                                                           'quality', 'heterozygosity'])
	paginator = Paginator(results, 50)
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)
		print "Error"

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/library_chromosome_filter.html', {"results": results,
	                                                                   "chromosome": chromosome,
	                                                                   "library": library,
	                                                                   "filter_urls": filter_urls,
	                                                                   "paginator": paginator,
	                                                                   "toolbar_max": toolbar_max,
	                                                                   "toolbar_min": toolbar_min,
	                                                                   "count": count})


# Displays the search page for a snp summary by library and chromosome level.
def chromosome_library_snp_summary(request):
	libraries = Library.objects.values('library_code').distinct().order_by('library_code')
	paginator = Paginator(libraries, 300)

	page = request.GET.get('page')
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)
	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/chromosome_library_snp_summary.html', {"results": results,
	                                                                        "paginator": paginator,
	                                                                        "toolbar_max": toolbar_max,
	                                                                        "toolbar_min": toolbar_min})


# Returns a chromosome level summary for an individual library.
def chromosome_library_snp_summary_filter(request):
	order_by = request.GET.get('order_by', 'chromosome__chromosome_name')
	library = request.GET.get('lib')
	results = SNP.objects.values('chromosome__chromosome_name',
	                             'library_id', 'library__library_code',
	                             'result__genome__organism__organismcode').filter(library__library_code=library).annotate(num_snps=Count('snp_id'),
	                                                                                                                      hetero=BooleanSum('heterozygosity'),
	                                                                                                                      indel=BooleanSum('snp_type__indel'),
	                                                                                                                      trans=BooleanSum('snp_type__transition'))
	organismcode = SNP.objects.values_list('result__genome__organism__organismcode', flat=True).filter(library__library_code=library)[0]
	library_size = get_chromosome_size(organismcode)
	result_list = results.order_by(order_by)
	paginator = Paginator(result_list, 50)
	page = request.GET.get('page')

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['chromosome__chromosome_name', 'library__librarysize', 'library_id', 'library__library_code',
	                                                           'num_snps', 'hetero', 'trans', 'indel'])

	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/chromosome_library_snp_summary_filter.html', {"results": results,
	                                                                               "library_size": library_size,
	                                                                               "filter_urls": filter_urls,
	                                                                               "paginator": paginator,
	                                                                               "toolbar_max": toolbar_max,
	                                                                               "toolbar_min": toolbar_min})


# Returns the snps found within a library and chromosome
def library_chromosome_snps(request):
	chromosomes = Chromosome.objects.values('chromosome_name').distinct().order_by('chromosome_name')
	page = request.GET.get('page')

	paginator = Paginator(chromosomes, 50)

	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)
	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/library_chromosome_snp.html', {"results": results,
	                                                                "paginator": paginator,
	                                                                "toolbar_max": toolbar_max,
	                                                                "toolbar_min": toolbar_min}, context_instance=RequestContext(request))


# todo change to access genes through the Effect table.
# Returns a full list of genes found within a specific library. Currently connects through the feature table.
def gene_list(request):
	order_by = request.GET.get('order_by', 'chromosome')
	library = request.GET.get('lib')
	results = Feature.objects.values('geneid', 'fmin', 'fmax', 'chromosome').filter(featuretype='gene',
	                                                                                genome__organism__library__library_code=library).order_by(order_by)
	count = len(results)

	paginator = Paginator(results, 200)
	page = request.GET.get('page')

	filter_urls = build_orderby_urls(request.get_full_path(), ['geneid', 'chromosome', 'fmin', 'fmax'])
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)
	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)

	return render_to_response('snpdb/gene_list.html', {"results": results,
	                                                   "library": library,
	                                                   "count": count,
	                                                   "filter_urls": filter_urls,
	                                                   "paginator": paginator,
	                                                   "toolbar_max": toolbar_max,
	                                                   "toolbar_min": toolbar_min})


# Displays the search page to compare two libraries for unique and similar snps.
def compare_two_libraries(request):
	lib_list = Library.objects.values('library_code', 'result__genome__organism__organismcode').distinct().order_by('result__genome__organism__organismcode')
	lib_genome = defaultdict(list)
	for each in lib_list:
		lib_code = str(each['library_code'])
		genome = str(each['result__genome__organism__organismcode'])
		if genome in lib_genome:
			cur_libs = lib_genome[genome]
			cur_libs.append(lib_code)
			lib_genome[genome] = cur_libs
		else:
			lib_codes = [lib_code]
			lib_genome[genome] = lib_codes
	page = request.GET.get('page')
	paginator = Paginator(tuple(lib_genome.items()), 500)

	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)
	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)
	return render_to_response('snpdb/compare_two_libraries.html', {"results": results,
	                                                               "paginator": paginator,
	                                                               "toolbar_max": toolbar_max,
	                                                               "toolbar_min": toolbar_min}, context_instance=RequestContext(request))


#todo flag snps with multiple alternates
# Lists identified snps (those found to be unique for a library) by their impact type. Passed from difference_two_libraries
def impact_snps(request):
	path = request.GET.get('path')
	library1 = request.GET.getlist('lib1')
	library2 = request.GET.getlist('lib2')
	lib = request.GET.get('library')

	impact = request.GET.get('impact')
	add_code = request.GET.get('samp1')
	neg_code = request.GET.get('samp2')
	order_by = request.GET.get('order_by', '0')
	print type(order_by)

	if add_code:
		cmd = """cat %s | /usr/local/Cellar/snpeff/3.6c/share/scripts/vcfEffOnePerLine.pl | java -jar /usr/local/Cellar/snpeff/3.6c/libexec/SnpSift.jar filter "( EFF[*].IMPACT = '%s' & !GEN[0].GT = '0/0')" | java -jar /usr/local/Cellar/snpeff/3.6c/libexec/SnpSift.jar extractFields - POS REF ALT CHROM EFF[*].GENE EFF[*].EFFECT QUAL"""
		library = library1
	elif neg_code:
		cmd = """cat %s | /usr/local/Cellar/snpeff/3.6c/share/scripts/vcfEffOnePerLine.pl | java -jar /usr/local/Cellar/snpeff/3.6c/libexec/SnpSift.jar filter "( EFF[*].IMPACT = '%s' & !GEN[1].GT = '0/0')" | java -jar /usr/local/Cellar/snpeff/3.6c/libexec/SnpSift.jar extractFields - POS REF ALT CHROM EFF[*].GENE EFF[*].EFFECT QUAL"""
		library = library2
	else:
		print "no sample passed"
		c ={"path": path,
	    "library1": library1,
	    "library2": library2}
		return render_to_response('snpdb/impact_snps_search.html', c, context_instance=RequestContext(request))

	snps_effect = subprocess.Popen(cmd % (path, impact), shell=True, stdout=subprocess.PIPE)

	snp_dict = defaultdict(dict)
	for line in snps_effect.stdout:
		snp = defaultdict(list)
		if line:
			if '#POS' not in line:
				entry = line.split('\t')
				ref = entry[1]
				alt = entry[2]
				pos = int(entry[0])
				qual = float(entry[6])
				gene = entry[4]
				chrom = entry[3].split('_')[0]
				impact = entry[5]
				effect = Feature.objects.filter(chromosome__startswith=chrom, featuretype='gene', geneid=gene).values('geneproduct')

				if pos in snp_dict:
					curr = snp_dict[pos]
					if ref not in curr['ref'] or alt not in curr['alt']:
						snp_dict[pos]['ref'].append(ref)
						snp_dict[pos]['alt'].append(alt)
					snp_dict[pos]['quality'].append(qual)
					snp_dict[pos]['gene'].append(gene)
					snp_dict[pos]['impact'].append(impact)
					snp_dict[pos]['effect'].append(str(effect[0]['geneproduct']))
				else:
					snp['chromosome'] = chrom
					snp['ref'] = [entry[1]]
					snp['alt'] = [entry[2]]
					snp['quality'] = [qual]
					snp['impact'] = [entry[5]]
					snp['gene'] = [gene]
					try:
						snp['effect'] = [str(effect[0]['geneproduct'])]
						entry.append(str(effect[0]['geneproduct']))
					except IndexError:
						pass
					snp_dict[pos] = snp
	if order_by == '0':
		sorted_snp = sorted(snp_dict.iteritems())
	else:
		sorted_snp = sorted(snp_dict.iteritems(), key=lambda (k, v): v[order_by])
	count = len(sorted_snp)
	paginator = Paginator(sorted_snp, 200)
	page = request.GET.get('page')

	# Calls utils method to append new filters or order_by to the current url
	filter_urls = build_orderby_urls(request.get_full_path(), ['0', 'chromosome', 'ref', 'alt', 'quality', 'impact', ])
	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)

	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)
	c ={"path": path,
	    "paginator": paginator,
	    "results": results,
	    "library": library,
	    "lib": lib,
	    "filter_urls": filter_urls,
	    "toolbar_max": toolbar_max,
	    "toolbar_min": toolbar_min,
	    "count": count, }
	return render_to_response('snpdb/impact_snps_search.html', c, context_instance=RequestContext(request))


# Compares multiple libraries by running bcftools isec. Returns the results and counts the number of snps by impact types.
def effects_by_vcf(request):
	library_1 = request.GET.getlist('check1')
	library_2 = request.GET.getlist('check2')

	#Captures vcf file location
	vcf1 = VCF_Files.objects.values_list('vcf_path', flat=True).filter(library__library_code__in=library_1)
	vcf2 = VCF_Files.objects.values_list('vcf_path', flat=True).filter(library__library_code__in=library_2).distinct()

	#Gets path of vcf files.
	direct = os.path.abspath(os.path.dirname(__file__))
	pro_dir = re.findall('(^.*)\/snpdb' , direct)[0]
	vcf_path = os.path.join(direct, 'vcf_files')



	#Collects all libraries to be compared.
	group_1_path = []
	for each in vcf1:
		vcf1_path = os.path.join(pro_dir, each)
		group_1_path.append(vcf1_path)
	group_2_path = []
	for each in vcf2:
		vcf2_path = os.path.join(pro_dir, each)
		group_2_path.append(vcf2_path)
	libraries = group_1_path + group_2_path
	libs = library_1 + library_2
	vcf_string = '_'.join(libs)

	#Determines the location of where analysis results will be stored.
	path = os.path.join(vcf_path, 'vcf_contrast_%s_%s' % (vcf_string, datetime.datetime.utcnow().strftime("%Y-%m-%d")))

	#Checks to see if analysis has already been completed. If the analysis files are not present, bcftools is called.
	if os.path.isdir(path):
		print "File already present"
		vcf_contrast = os.path.join(path, 'vcf_contrast.vcf')
		pass
	else:
		os.mkdir(path)
		#Checks to see if files have been zipped and indexed. Bcftools requires indexed vcf files.
		for fname in libraries:
			if os.path.isfile(fname):
				# zips and indexes vcf-files for vcf-merge
				try:
					subprocess.check_call(['bgzip', fname])
					subprocess.check_call(['tabix', '-p', 'vcf', '%s.gz' % fname])
				except IOError:
					pass
			elif os.path.isfile('%s.gz' % fname):
				print "files already zipped"

		#Creates a string of all zipped file
		zip_vcf = ''
		for each in libraries:
			zips = str(each) + '.gz'
			zip_vcf = zip_vcf + ' ' + zips

		#Runs the vcf-merge command.
		merge_file = os.path.join(path, 'merge.vcf')
		p = subprocess.Popen(["""/usr/local/bin/bcftools merge --force-samples %s > %s""" % (zip_vcf, merge_file)], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
		out, err = p.communicate()
		print "files merged"

		#collects all file ids.
		add_code = []
		neg_code = []

		h = subprocess.Popen(["""/usr/bin/grep ^#CHROM %s""" % merge_file], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
		header, error = h.communicate()
		replace = '\t'.join(libs).strip('\n')
		samp = re.findall(r'^#CHROM\t\w*\t*\w*\t\w*\t\w*\t\w*\t\w*\t\w*\t\w*\t(.*)', header)[0]
		replace_str = re.sub(samp, replace, header).strip('\n')
		sed = subprocess.call(["""/usr/bin/perl -pi -e 's/^#CHROM.*/%s/g;' %s""" % (replace_str, merge_file)], shell=True)
		print "sed is complete"


		for lib in library_1:
			add_code.append(lib)
		for lib in library_2:
			neg_code.append(lib)
		add = '+' + ','.join(add_code)
		neg = '-' + ','.join(neg_code)

		proj_dir = os.path.abspath(os.path.join(direct, os.pardir))
		script_dir = os.path.join(proj_dir, 'scripts')

		#splits the multiallelic entries in the merge file
		split_command = os.path.join(script_dir, 'split_multi_add_WT.py')
		replace_file = os.path.join(path, 'merge_split.vcf')
		subprocess.call(["""python %s %s %s""" % (split_command, merge_file, replace_file)], shell=True)
		print "multi-allelic sites split"

		#removes all blank entries and replaces them with WT stats.
		# replace_command = os.path.join(script_dir, 'add_standard_WT.py')
		# replace_file = os.path.join(path, 'merge_split_replace.vcf')
		# subprocess.call(["""/usr/bin/python %s %s %s""" % (replace_command, split_file, replace_file)], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
		# print "replaced empty entries with WT"

		#zips replace file for vcf-contrast
		try:
			subprocess.check_call(['bgzip', replace_file])
			subprocess.check_call(['tabix', '-p', 'vcf', '%s.gz' % replace_file])
		except IOError:
			pass

		print "calling vcf-contrast"
		zip_replace = replace_file + '.gz'
		vcf_contrast = os.path.join(path, 'vcf_contrast.vcf')
		subprocess.call(["""/usr/local/bin/vcf-contrast -n %s %s %s > %s""" % (add, neg, zip_replace, vcf_contrast)], shell=True)

	#Opens the returned vcf-contrast file and counts the data.
	print "opening vcf-contrast"
	lib1_effect = []
	lib2_effect = []
	lib1_total = 0
	lib2_total = 0

	vcf_reader = vcf.Reader(open ('%s' % vcf_contrast, 'r'))
	date = datetime.datetime.utcnow().strftime("%Y-%m-%d").replace('-', '')
	source = 'source_' + date + '.1'
	cmd = vcf_reader.metadata[source][0]
	add = re.findall('\+(.*) \-', cmd)[0].split(',')
	neg = re.findall('\+.* \-(.*) ', cmd)[0].split(',')

	lib1_high_effects = defaultdict(int)
	lib1_moderate_effects = defaultdict(int)
	lib1_modifier_effects = defaultdict(int)
	lib1_low_effects = defaultdict(int)

	lib2_high_effects = defaultdict(int)
	lib2_moderate_effects = defaultdict(int)
	lib2_modifier_effects = defaultdict(int)
	lib2_low_effects = defaultdict(int)

	lib1_total_counts = [0, 0, 0, 0, 0, '']
	lib2_total_counts = [0, 0, 0, 0, 0, '']
	for record in vcf_reader:
		lib1 = False
		lib2 = False
		#Keeps track of what effect type each snp has. [high, moderate, low, modifier]
		lib1_impact_counts = [0, 0, 0, 0]
		lib2_impact_counts = [0, 0, 0, 0]
		#Places each type of impact into dictionary of the effect. SNPs with multiple impacts will have all impacts accounted for in the impact total.
		# i.e, SNPs with Downstream and Upstream effects will results in an addition to both impact counts.
		# for x in effects:
		for sample in record.samples:
			samp_id = sample.sample
			gt = sample['GT']
			if gt != '0/0':
				if samp_id in add:
					lib1 = True
				if samp_id in neg:
					lib2 = True
		effects = record.INFO['EFF']
		for x in effects:
			impact = x.split('(')[0]
			effect_list = x.split('(')[1]
			effect = effect_list.split('|')
			if effect[0] == "HIGH":
				if lib1 is True:
					lib1_impact_counts[0] += 1
					lib1_high_effects[impact] += 1
				if lib2 is True:
					lib2_impact_counts[0] += 1
					lib2_high_effects[impact] += 1
			elif effect[0] == "MODERATE":
				if lib1 is True:
					lib1_impact_counts[1] += 1
					lib1_moderate_effects[impact] += 1
				if lib2 is True:
					lib2_impact_counts[1] += 1
					lib2_moderate_effects[impact] += 1
			elif effect[0] == "LOW":
				if lib1 is True:
					lib1_impact_counts[2] += 1
					lib1_low_effects[impact] += 1
				if lib2 is True:
					lib2_impact_counts[2] += 1
					lib2_low_effects[impact] += 1
			elif effect[0] == "MODIFIER":
				if lib1 is True:
					lib1_impact_counts[3] += 1
					lib1_modifier_effects[impact] += 1
				if lib2 is True:
					lib2_impact_counts[3] += 1
					lib2_modifier_effects[impact] += 1

		#Counts the number of snps effected by each impact type. Snp is only counted once for each impact, i.e. if SNP has two modifying impacts, it is only counted once.
		if sum(lib1_impact_counts) > 0:
			lib1_total_counts[5] = add
			lib1_total_counts[4] += 1

		if lib1_impact_counts[0] > 0:
			lib1_total_counts[0] += 1
		if lib1_impact_counts[1] > 0:
			lib1_total_counts[1] += 1
		if lib1_impact_counts[2] > 0:
			lib1_total_counts[2] += 1
		if lib1_impact_counts[3] > 0:
			lib1_total_counts[3] += 1

		if sum(lib2_impact_counts) > 0:
			lib2_total_counts[5] = neg
			lib2_total_counts[4] += 1
		if lib2_impact_counts[0] > 0:
			lib2_total_counts[0] += 1
		if lib2_impact_counts[1] > 0:
			lib2_total_counts[1] += 1
		if lib2_impact_counts[2] > 0:
			lib2_total_counts[2] += 1
		if lib2_impact_counts[3] > 0:
			lib2_total_counts[3] += 1
	lib1_tuple = (dict(lib1_high_effects), dict(lib1_moderate_effects), dict(lib1_modifier_effects), dict(lib1_low_effects), lib1_total_counts)
	lib1_effect.append(lib1_tuple)
	lib1_total += lib1_total_counts[4]

	lib2_tuple = (dict(lib2_high_effects), dict(lib2_moderate_effects), dict(lib2_modifier_effects), dict(lib2_low_effects), lib2_total_counts)
	lib2_effect.append(lib2_tuple)
	lib2_total += lib2_total_counts[4]

	return render_to_response('snpdb/test2.html', {"library1": library_1,
	                                               "library2": library_2,
	                                               "lib1_high_effects": dict(lib1_high_effects),
	                                               "lib1_moderate_effects": dict(lib1_moderate_effects),
	                                               "lib1_modifier_effects": dict(lib1_modifier_effects),
	                                               "lib1_low_effects": dict(lib1_low_effects),
	                                               "lib1_total_counts": lib1_total_counts,
	                                               "lib1_total": lib1_total,
	                                               "lib1_effect": lib1_effect,
	                                               "lib2_high_effects": dict(lib2_high_effects),
	                                               "lib2_moderate_effects": dict(lib2_moderate_effects),
	                                               "lib2_modifier_effects": dict(lib2_modifier_effects),
	                                               "lib2_low_effects": dict(lib2_low_effects),
	                                               "lib2_total_counts": lib2_total_counts,
	                                               "lib2_total": lib2_total,
	                                               "lib2_effect": lib2_effect,
	                                               "add_code": add,
	                                               "neg_code": neg,
	                                               "path": vcf_contrast})


def get_chromosome_size(organismcode):
	genome_size = Chromosome.objects.filter(genome_name_id=organismcode).aggregate(genome_size=Sum('size'))
	return genome_size


# def impact_snps2(request):
# 	order_by = request.GET.get('order_by', 'snp_id').encode("ascii")
# 	impact = request.GET.get('impact')
# 	library1 = request.GET.get('lib1')
# 	library2 = request.GET.get('lib2')
#
# 	snp2 = SNP.objects.filter(library__library_code=library2).values_list('snp_position',
# 	                                                                      'chromosome__chromosome_name')
# 	snp1 = SNP.objects.filter(library__library_code=library1).values_list('snp_position',
# 	                                                                      'chromosome__chromosome_name')
# 	difference = set(snp1).difference(set(snp2))
# 	snps = []
# 	for x in difference:
# 		keys = ['snp_position', 'chromosome__chromosome_name']
# 		snps.append(dict(zip(keys, x)))
#
# 	snp_ids = []
# 	for each in snps:
# 		ids = SNP.objects.values_list('snp_id', flat=True).filter(library__library_code=library1, snp_position=each['snp_position'],
# 		                                                          chromosome__chromosome_name=each['chromosome__chromosome_name'])
# 		snp_ids.append(ids)
# 	snp_id = [x for sublist in snp_ids for x in sublist]
#
# 	snps = SNP.objects.filter(snp_id__in=snp_id, effect__effect_string=impact).values_list('snp_id', 'effect__effect_group')
# 	groups = []
# 	for y in snps:
# 		keys = ['snp_id', 'effect__effect_group']
# 		groups.append(dict(zip(keys, y)))
# 	genes = []
# 	for each in groups:
# 		genes.append(SNP.objects.filter(**each).values('library', 'library__library_code', 'snp_id',
# 		                                               'snp_position', 'ref_base', 'alt_base',
# 		                                               'heterozygosity', 'quality', 'result_id',
# 		                                               'chromosome__chromosome_name', 'effect__effect_string',
# 		                                               'effect__effect_class', 'effect__effect'))
# 	group = [x for sublist in genes for x in sublist]
# 	sorted_gene_dict = genes_from_effect(group, library1, order_by)
# 	gene_dict = [x for sublist in sorted_gene_dict for x in sublist]
# 	print len(filter(None, gene_dict))
# 	# new_dict = {k:v for k,v in gene_dict.items() if v}
# 	count = len(gene_dict)
# 	paginator = Paginator(gene_dict, 200)
# 	page = request.GET.get('page')
#
# 	# Calls utils method to append new filters or order_by to the current url
# 	filter_urls = build_orderby_urls(request.get_full_path(), ['library', 'library__library_code', 'snp_id',
# 	                                                           'snp_position', 'ref_base', 'alt_base',
# 	                                                           'heterozygosity', 'quality',
# 	                                                           'chromosome__chromosome_name',
# 	                                                           'effect__effect_string'])
# 	try:
# 		results = paginator.page(page)
# 	except PageNotAnInteger:
# 		results = paginator.page(1)
# 	except EmptyPage:
# 		results = paginator.page(paginator.num_pages)
#
# 	toolbar_max = min(results.number + 4, paginator.num_pages)
# 	toolbar_min = max(results.number - 4, 0)
# 	return render_to_response('snpdb/impact_snps_search.html', {"results": results,
# 	                                                            "toolbar_max": toolbar_max,
# 	                                                            "toolbar_min": toolbar_min,
# 	                                                            "count": count,
# 	                                                            "filter_urls": filter_urls,
# 	                                                            })


def chrom_region(request):
	lib_list = Library.objects.values('library_code', 'result__genome__organism__organismcode').order_by('result__genome__organism__organismcode')
	lib_genome = defaultdict(list)
	for each in lib_list:
		lib_code = str(each['library_code'])
		genome = str(each['result__genome__organism__organismcode'])
		if genome in lib_genome:
			cur_libs = lib_genome[genome]
			cur_libs.append(lib_code)
			lib_genome[genome] = cur_libs
		else:
			lib_codes = [lib_code]
			lib_genome[genome] = lib_codes
	page = request.GET.get('page')
	paginator = Paginator(tuple(lib_genome.items()), 500)

	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)
	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)
	return render_to_response('snpdb/chrom_region.html', {"results": results,
	                                                      "paginator": paginator,
	                                                      "toolbar_max": toolbar_max,
	                                                      "toolbar_min": toolbar_min}, context_instance=RequestContext(request))


def chrom_region_search(request):
	library = request.GET.get('check1')
	chromosome = SNP.objects.values_list('chromosome__chromosome_name').filter(library__library_code=library).distinct().order_by('chromosome__chromosome_name')
	return render_to_response('snpdb/chrom_region.html', {"chromosome": chromosome,
	                                                      "library": library,
	                                                      })


def chrom_region_filter(request):
	library = request.GET.get('library')
	start = request.GET.get('from')
	stop = request.GET.get('to')
	chrom = request.GET.get('chrom')

	genes = SNP.objects.filter(effect__effect_id=6).values('snp_position', 'effect__effect_string',
	                                                       'effect__effect_group').filter(chromosome__chromosome_name=chrom, library__library_code=library,
	                                                                                      snp_position__range=(start, stop)).distinct()
	results = defaultdict(int)
	for each in genes:
		effect_group = each['effect__effect_group']
		impact = SNP.objects.filter(effect__effect_id=1,
		                            effect__effect_group=effect_group,).values('snp_position', 'ref_base', 'alt_base', 'heterozygosity',
		                                                                       'quality', 'effect__effect_string', 'effect__effect_group').filter(
			snp_position=each['snp_position'],
			chromosome__chromosome_name=chrom,
			library__library_code=library,).distinct()
		for x in impact:
			x['gene'] = each['effect__effect_string']
			# print "x: ", x
			if each['snp_position'] in results:
				# print "current: ", each
				results[each['snp_position']].append(x)
			else:
				results[each['snp_position']] = [x,]
			# print "results: ", results[each['snp_position']]

	page = request.GET.get('page')
	paginator = Paginator(results.items(), 50)

	try:
		results = paginator.page(page)
	except PageNotAnInteger:
		results = paginator.page(1)
	except EmptyPage:
		results = paginator.page(paginator.num_pages)
	toolbar_max = min(results.number + 4, paginator.num_pages)
	toolbar_min = max(results.number - 4, 0)
	print toolbar_min, toolbar_max
	return render_to_response('snpdb/chrom_region_filter.html', {"chromosome": chrom,
	                                                             "library": library,
	                                                             "results": results,
	                                                             "toolbar_max": toolbar_max,
	                                                             "toolbar_min": toolbar_min}, context_instance=RequestContext(request))


# Dumps a queryset into a csv file.
def dump(qs, outfile_path):
	writer = csv.writer(open(outfile_path, 'w'))
	keys = []
	values = []
	for obj in qs:
		for key, val in obj.items():
			values.append(val)
			if key not in keys:
				keys.append(key)
	value = [tuple(values[i:i+2]) for i in range(0, len(values), 2)]
	writer.writerow(keys)
	for each in value:
		writer.writerow(each)


# Reads a csv file and converts the data back into a dictionary.
def read(filename):
	impact_dict = {}
	for each in csv.reader(open(filename)):
		for x in each:
			x = x.split(',')
			key = x[0].replace('(', '').replace("'", '')
			value = x[1].strip(')').replace("'", '')
			try:
				impact_dict[key] = int(value)
			except ValueError:
				impact_dict[key] = value
				pass
	return impact_dict


# unused code.
#----------------------------------------------------------------------------------------------------------------------

# Commands to save the snpdb dashboard pie-charts. Should be run after each vcf import.
def save_snp_dashboard_files(chart_path, image_path):
	# Google Chart Images
	lib_labels = []
	lib_legend = []
	org_labels = []
	org_legend = []
	impact_labels = []
	high_labels = []
	low_labels = []
	moderate_labels = []
	modifier_labels = []
	high_keys = []
	low_keys = []
	moderate_keys = []
	modifier_keys = []
	impact_keys = []
	high_values = []
	low_values = []
	moderate_values = []
	modifier_values = []
	impact_values = []

	high = Effect.objects.filter(effect_id=1, effect_string="HIGH").values("effect_class").annotate(Count('snp'))
	dump(high, chart_path % 'high')
	for obj in high.iterator():
		for key, val in obj.items():
			high_values.append(val)
			if val not in high_keys and not isinstance(val, int):
				high_keys.append(val)
	high_value = [tuple(high_values[i:i+2]) for i in range(0, len(high_values), 2)]
	high_snp_total = sum(i[1] for i in high_value)
	for x in high_value:
		percentage = float(x[1])/float(high_snp_total)*100
		high_labels.append(round(percentage, 2))
	snps_by_high_impact = Pie(high_labels).label(*high_labels).legend(*high_keys).color("919dab", "D2E3F7",
	                                                                                    "658CB9", "88BBF7",
	                                                                                    "666E78").size(450, 200)
	snps_by_high_impact.image().save(image_path % 'high', 'png')
	print "high files saved"

	impact = Effect.objects.filter(effect_id=1).values("effect_string").annotate(Count('snp'))
	dump(impact, chart_path % 'impact')
	for obj in impact.iterator():
		for key, val in obj.items():
			print key, val
			impact_values.append(val)
			if val not in impact_keys and not isinstance(val, int):
				impact_keys.append(val)
	impact_value = [tuple(impact_values[i:i+2]) for i in range(0, len(impact_values), 2)]
	impact_snp_total = sum(i[1] for i in impact_value)
	for x in impact_value:
		percentage = float(x[1])/float(impact_snp_total)*100
		impact_labels.append(round(percentage,2))
	snps_by_impact = Pie(impact_labels).label(*impact_labels).legend(*impact_keys).color("919dab", "D2E3F7",
	                                                                                     "658CB9", "88BBF7",
	                                                                                     "666E78").size(450, 200)
	snps_by_impact.image().save(image_path % 'impact', 'png')
	print "impact files saved"

	low = Effect.objects.filter(effect_id=1, effect_string="LOW").values("effect_class").annotate(Count('snp'))
	dump(low, chart_path % 'low')
	for obj in low.iterator():
		for key, val in obj.items():
			low_values.append(val)
			if val not in low_keys and not isinstance(val, int):
				low_keys.append(val)
	low_value = [tuple(low_values[i:i+2]) for i in range(0, len(low_values), 2)]
	low_snp_total = sum(i[1] for i in low_value)
	for x in low_value:
		percentage = float(x[1])/float(low_snp_total)*100
		low_labels.append(round(percentage, 2))
	snps_by_low = Pie(low_labels).label(*low_labels).legend(*low_keys).color("919dab", "D2E3F7",
	                                                                         "658CB9", "88BBF7",
	                                                                         "666E78").size(450, 200)
	snps_by_low.image().save(image_path % 'low', 'png')
	print "low files saved"

	moderate = Effect.objects.filter(effect_id=1, effect_string="MODERATE").values("effect_class").annotate(Count('snp'))
	dump(moderate, chart_path % 'moderate')
	for obj in moderate.iterator():
		for key, val in obj.items():
			moderate_values.append(val)
			if val not in moderate_keys and not isinstance(val, int):
				moderate_keys.append(val)
	moderate_value = [tuple(moderate_values[i:i+2]) for i in range(0, len(moderate_values), 2)]
	moderate_snp_total = sum(i[1] for i in moderate_value)
	for x in moderate_value:
		percentage = float(x[1])/float(moderate_snp_total)*100
		moderate_labels.append(round(percentage, 2))
	snps_by_moderate = Pie(moderate_labels).label(*moderate_labels).legend(*moderate_keys).color("919dab", "D2E3F7",
	                                                                                             "658CB9", "88BBF7",
	                                                                                             "666E78").size(550, 200)
	snps_by_moderate.image().save(image_path % 'moderate', 'png')
	print "moderate files saved"

	modifier = Effect.objects.filter(effect_id=1, effect_string="MODIFIER").values("effect_class").annotate(Count('snp'))
	dump(modifier, chart_path % 'modifier')
	for obj in modifier.iterator():
		for key, val in obj.items():
			modifier_values.append(val)
			if val not in modifier_keys and not isinstance(val, int):
				modifier_keys.append(val)
	modifier_value = [tuple(modifier_values[i:i+2]) for i in range(0, len(modifier_values), 2)]
	modifier_snp_total = sum(i[1] for i in modifier_value)
	for x in modifier_value:
		percentage = float(x[1])/float(modifier_snp_total)*100
		modifier_labels.append(round(percentage, 2))
	snps_by_modifier = Pie(modifier_labels).label(*modifier_labels).legend(*modifier_keys).color("919dab", "D2E3F7",
	                                                                                             "658CB9", "88BBF7",
	                                                                                             "666E78").size(450, 200)
	snps_by_modifier.image().save(image_path % 'modifier', 'png')
	print "modifier files saved"

	lib_count = SNP.objects.values("library__library_code").distinct().annotate(Count('snp_id'))
	lib_snps = []
	lib_snp_total = 0
	for each in lib_count.iterator():
		lib_snps.append(each['snp_id__count'])
		lib_snp_total += each['snp_id__count']
	for x in lib_snps:
		percentage = float(x)/float(lib_snp_total)*100
		lib_labels.append(round(percentage, 2))
	lib_legend.append(each['library__library_code'])
	snps_by_library = Pie([lib_labels]).label(*lib_labels).legend(*lib_legend).color("919dab", "D2E3F7",
	                                                                                 "658CB9", "88BBF7",
	                                                                                 "666E78").size(450, 200)
	snps_by_library.image()
	snps_by_library.image().save(image_path % 'library', 'png')
	print "saved snps_by_library"

	org_count = SNP.objects.values("library__organism__organismcode").distinct().annotate(Count('snp_id'))
	org_snps = []
	org_snp_total = 0
	for each in org_count.iterator():
		org_snps.append(each['snp_id__count'])
		org_snp_total += each['snp_id__count']
	for x in org_snps:
		percentage = float(x)/float(org_snp_total)*100
		org_labels.append(round(percentage, 2))
	org_legend.append(each['library__organism__organismcode'])
	snps_by_organism = Pie(org_labels).label(*org_labels).legend(*org_legend).color("919dab", "D2E3F7",
	                                                                                "658CB9", "88BBF7",
	                                                                                "666E78").size(450, 200)
	snps_by_organism.image().save(image_path % 'organism', 'png')
	print "saved snps_by_organism"



# View to search a list of genes for snps.
# def multi_gene_snps(request):
# 	return render_to_response('snpdb/multi_gene_snps.html',)


# Returns snps found within a list of genes. Does not consider library.
# def multi_gene_snps_filter(request):
# 	order_by = request.GET.get('order_by', 'effect__effect_string')
# 	gene = request.GET.get('s')
# 	genes = gene.split()
#
# 	result_list = SNP.objects.filter(effect__effect_id=6, effect__effect_string__in=genes,
# 	                                 effect__effect_class__endswith='SYNONYMOUS_CODING'.decode('utf-8')).values('library', 'library__library_code', 'snp_id',
# 	                                                                                                            'snp_position', 'ref_base', 'alt_base',
# 	                                                                                                            'heterozygosity', 'quality',
# 	                                                                                                            'chromosome__chromosome_name', 'effect__effect_string',
# 	                                                                                                            'effect__effect_class', 'effect__effect', 'result_id').order_by(order_by)
#
# 	count = result_list.count()
# 	paginator = Paginator(result_list, 50)
# 	page = request.GET.get('page')
# 	filter_urls = build_orderby_urls(request.get_full_path(), ['library', 'library__library_code', 'snp_id',
# 	                                                           'snp_position', 'ref_base', 'alt_base',
# 	                                                           'heterozygosity', 'quality',
# 	                                                           'chromosome__chromosome_name', 'effect__effect_string',
# 	                                                           'effect__effect_class', 'effect__effect', 'result_id'])
# 	try:
# 		results = paginator.page(page)
# 	except PageNotAnInteger:
# 		results = paginator.page(1)
# 	except EmptyPage:
# 		results = paginator.page(paginator.num_pages)
#
# 	toolbar_max = min(results.number + 4, paginator.num_pages)
# 	toolbar_min = max(results.number - 4, 0)
#
# 	return render_to_response('snpdb/multi_gene_snps_filter.html', {"results": results,
# 	                                                                "filter_urls": filter_urls,
# 	                                                                "paginator": paginator,
# 	                                                                "toolbar_max": toolbar_max,
# 	                                                                "toolbar_min": toolbar_min,
# 	                                                                "genes": genes,
# 	                                                                "count": count})


# def multi_gene_library_snps(request):
# 	libraries = Library.objects.values('library_code').order_by('library_code')
# 	return render_to_response('snpdb/multi_gene_snps_library.html', {"libraries": libraries, })


# def multi_gene_library_snps_filter(request):
# 	order_by = request.GET.get('order_by', 'effect__effect_string')
# 	gene = request.GET.get('s')
# 	libraries = request.GET.getlist('check')
# 	genes = gene.split()
# 	result_list = SNP.objects.filter(effect__effect_id=6, effect__effect_string__in=genes,
# 	                                 library__library_code__in=libraries).values('library', 'library__library_code', 'snp_id',
# 	                                                                             'snp_position', 'ref_base', 'alt_base',
# 	                                                                             'heterozygosity', 'quality',
# 	                                                                             'chromosome__chromosome_name', 'effect__effect_string',
# 	                                                                             'effect__effect_class', 'effect__effect', 'result_id').order_by(order_by)
#
# 	paginator = Paginator(result_list, 50)
# 	page = request.GET.get('page')
# 	filter_urls = build_orderby_urls(request.get_full_path(), ['library', 'library__library_code', 'snp_id',
# 	                                                           'snp_position', 'ref_base', 'alt_base',
# 	                                                           'heterozygosity', 'quality',
# 	                                                           'chromosome__chromosome_name', 'effect__effect_string',
# 	                                                           'effect__effect_class', 'effect__effect', 'result_id'])
# 	try:
# 		results = paginator.page(page)
# 	except PageNotAnInteger:
# 		results = paginator.page(1)
# 	except EmptyPage:
# 		results = paginator.page(paginator.num_pages)
#
# 	toolbar_max = min(results.number + 4, paginator.num_pages)
# 	toolbar_min = max(results.number - 4, 0)
#
# 	return render_to_response('snpdb/multi_gene_snps_library_filter.html', {"results": results,
# 	                                                                        "filter_urls": filter_urls,
# 	                                                                        "paginator": paginator,
# 	                                                                        "toolbar_max": toolbar_max,
# 	                                                                        "toolbar_min": toolbar_min,
# 	                                                                        "genes": genes})


#Identifies SNPs present in one library that are not present in the second. Uses python set differences.
# def diff_libraries2(request):
# 	library1 = request.GET.get('lib1')
# 	library2 = request.GET.get('lib2')
#
# 	snp2 = SNP.objects.filter(library__library_code=library2).values_list('snp_position',
# 	                                                                      'chromosome__chromosome_name')
# 	snp1 = SNP.objects.filter(library__library_code=library1).values_list('snp_position',
# 	                                                                      'chromosome__chromosome_name')
# 	difference = set(snp1).difference(set(snp2))
# 	opp_diff = set(snp2).difference(set(snp1))
# 	snps = []
# 	for x in difference:
# 		keys = ['snp_position', 'chromosome__chromosome_name']
# 		snps.append(dict(zip(keys, x)))
#
# 	snp_ids = []
# 	for each in snps:
# 		ids = SNP.objects.values_list('snp_id', flat=True).filter(library__library_code=library1, snp_position=each['snp_position'],
# 		                                                          chromosome__chromosome_name=each['chromosome__chromosome_name'])
# 		snp_ids.append(ids)
# 	snp_id = [x for sublist in snp_ids for x in sublist]
# 	opp_snps = []
# 	for x in opp_diff:
# 		keys = ['snp_position', 'chromosome__chromosome_name']
# 		opp_snps.append(dict(zip(keys, x)))
#
# 	opp_snp_ids = []
# 	for each in opp_snps:
# 		ids = SNP.objects.values_list('snp_id', flat=True).filter(library__library_code=library1, snp_position=each['snp_position'],
# 		                                                          chromosome__chromosome_name=each['chromosome__chromosome_name'])
# 		opp_snp_ids.append(ids)
# 	opp_snp_id = [x for sublist in opp_snp_ids for x in sublist]
#
# 	print "got result"
#
# 	snp_impact = Effect.objects.filter(snp__in=snp_id, effect=1).values('effect_string').annotate(snp_count=Count('snp')).order_by('effect_string')
# 	effects = Effect.objects.filter(snp__in=snp_id, effect=1).values('effect', 'effect_class',
# 	                                                                 'effect_string').annotate(effect_count=Count('snp')).order_by('effect_class')
# 	opp_snp_impact = Effect.objects.filter(snp__in=opp_snp_id, effect=1).values('effect_string').annotate(snp_count=Count('snp')).order_by('effect_string')
# 	opp_effects = Effect.objects.filter(snp__in=opp_snp_id, effect=1).values('effect', 'effect_class',
# 	                                                                         'effect_string').annotate(effect_count=Count('snp')).order_by('effect_class')
# 	print "got modifier"
# 	return render_to_response('snpdb/db_impact_snps.html', {"effects": effects,
# 	                                                        "opp_effects": opp_effects,
# 	                                                        "snp_impact": snp_impact,
# 	                                                        "opp_snp_impact": opp_snp_impact,
# 	                                                        "library1": library1,
# 	                                                        "library2": library2})


# Uses bcftools and SnpSift to identify snps that are unique between multiple libraries. This information is organized by impact type.
# def difference_two_libraries(request):
# 	library1 = request.GET.get('lib1')
# 	library2 = request.GET.get('lib2')
# 	vcf_path = os.path.join(dir, 'vcf_files')
# 	path = os.path.join(vcf_path, 'bcftools_isec_snpEff_%s_%s_%s' % (library1, library2, datetime.datetime.utcnow().strftime("%Y-%m-%d")))
# 	if os.path.isdir(path):
# 		print "file already present"
# 		pass
# 	else:
# 		print "initial analysis being ran"
# 		test = subprocess.check_call(["""bcftools isec %s/%s_gatk.snpEff.vcf.gz %s/%s_gatk.snpEff.vcf.gz -p %s""" % (vcf_path, library1, vcf_path, library2, path)],
# 		                             shell=True)
# 	cmd = """cat %s/0000.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 	modifier = subprocess.Popen(cmd % (path, "MODIFIER"), shell=True, stdout=subprocess.PIPE)
# 	moderate = subprocess.Popen(cmd % (path, "MODERATE"), shell=True, stdout=subprocess.PIPE)
# 	high = subprocess.Popen(cmd % (path, "HIGH"), shell=True, stdout=subprocess.PIPE)
# 	low = subprocess.Popen(cmd % (path, "LOW"), shell=True, stdout=subprocess.PIPE)
# 	total = subprocess.Popen("""grep ^[^#] %s/0000.vcf | wc -l""" % path, shell=True, stdout=subprocess.PIPE)
# 	counts = [high.communicate()[0].strip(), moderate.communicate()[0].strip(), low.communicate()[0].strip(), modifier.communicate()[0].strip(), total.communicate()[0].strip()]
#
# 	low_counts = defaultdict(int)
# 	high_counts = defaultdict(int)
# 	moderate_counts = defaultdict(int)
# 	modifier_counts = defaultdict(int)
# 	for each in low_effects:
# 		count_effect_cmd = """cat %s/0000.vcf| cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 		count_effect = subprocess.Popen(count_effect_cmd % (path, each), shell=True, stdout=subprocess.PIPE)
# 		count = count_effect.communicate()[0]
# 		low_counts[each] = count.strip()
# 	for each in moderate_effects:
# 		count_effect_cmd = """cat %s/0000.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 		count_effect = subprocess.Popen(count_effect_cmd % (path, each), shell=True, stdout=subprocess.PIPE)
# 		count = count_effect.communicate()[0]
# 		moderate_counts[each] = count.strip()
# 	for each in modifier_effects:
# 		count_effect_cmd = """cat %s/0000.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 		count_effect = subprocess.Popen(count_effect_cmd % (path, each), shell=True, stdout=subprocess.PIPE)
# 		count = count_effect.communicate()[0]
# 		modifier_counts[each] = count.strip()
# 	for each in high_effects:
# 		count_effect_cmd = """cat %s/0000.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 		count_effect = subprocess.Popen(count_effect_cmd % (path, each), shell=True, stdout=subprocess.PIPE)
# 		count = count_effect.communicate()[0]
# 		high_counts[each] = count.strip()
#
#
# 	cmd2 = """cat %s/0001.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 	modifier2 = subprocess.Popen(cmd2 % (path, "MODIFIER"), shell=True, stdout=subprocess.PIPE)
# 	moderate2 = subprocess.Popen(cmd2 % (path, "MODERATE"), shell=True, stdout=subprocess.PIPE)
# 	high2 = subprocess.Popen(cmd2 % (path, "HIGH"), shell=True, stdout=subprocess.PIPE)
# 	low2 = subprocess.Popen(cmd2 % (path, "LOW"), shell=True, stdout=subprocess.PIPE)
# 	total2 = subprocess.Popen("""grep ^[^#] %s/0001.vcf | wc -l""" % path, shell=True, stdout=subprocess.PIPE)
# 	counts2 = [high2.communicate()[0].strip(), moderate2.communicate()[0].strip(), low2.communicate()[0].strip(), modifier2.communicate()[0].strip(), total2.communicate()[0].strip()]
#
#
# 	low_counts2 = defaultdict(int)
# 	high_counts2 = defaultdict(int)
# 	moderate_counts2 = defaultdict(int)
# 	modifier_counts2 = defaultdict(int)
# 	for each in low_effects:
# 		count_effect_cmd = """cat %s/0001.vcf| cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 		count_effect = subprocess.Popen(count_effect_cmd % (path, each), shell=True, stdout=subprocess.PIPE)
# 		count = count_effect.communicate()[0]
# 		low_counts2[each] = count.strip()
# 	for each in moderate_effects:
# 		count_effect_cmd = """cat %s/0001.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 		count_effect = subprocess.Popen(count_effect_cmd % (path, each), shell=True, stdout=subprocess.PIPE)
# 		count = count_effect.communicate()[0]
# 		moderate_counts2[each] = count.strip()
# 	for each in modifier_effects:
# 		count_effect_cmd = """cat %s/0001.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 		count_effect = subprocess.Popen(count_effect_cmd % (path, each), shell=True, stdout=subprocess.PIPE)
# 		count = count_effect.communicate()[0]
# 		modifier_counts2[each] = count.strip()
# 	for each in high_effects:
# 		count_effect_cmd = """cat %s/0001.vcf | cut -f 8 | tr ";" "\n" | grep ^EFF= | cut -f 2 -d = | tr "," "\n" | grep %s | wc -l"""
# 		count_effect = subprocess.Popen(count_effect_cmd % (path, each), shell=True, stdout=subprocess.PIPE)
# 		count = count_effect.communicate()[0]
# 		high_counts2[each] = count.strip()
# 	return render_to_response('snpdb/impact_snps.html', {"counts": counts,
# 	                                                     "counts2": counts2,
# 	                                                     "low_counts": dict(low_counts),
# 	                                                     "high_counts": dict(high_counts),
# 	                                                     "moderate_counts": dict(moderate_counts),
# 	                                                     "modifier_counts": dict(modifier_counts),
# 	                                                     "low_counts2": dict(low_counts2),
# 	                                                     "high_counts2": dict(high_counts2),
# 	                                                     "moderate_counts2": dict(moderate_counts2),
# 	                                                     "modifier_counts2": dict(modifier_counts2),
# 	                                                     "library1": library1,
# 	                                                     "library2": library2,
# 	                                                     "path": path}, context_instance=RequestContext(request))
#


# Compares multiple libraries by running bcftools isec. Returns the results and counts the number of snps by impact types.
# def effects_by_vcf(request):
# 	library_1 = request.GET.getlist('check1')
# 	library_2 = request.GET.getlist('check2')
#
# 	#Captures vcf file location
# 	vcf1 = VCF_Files.objects.values_list('vcf_path', flat=True).filter(library__library_code__in=library_1)
# 	vcf2 = VCF_Files.objects.values_list('vcf_path', flat=True).filter(library__library_code__in=library_2).distinct()
#
# 	#Gets path of vcf files.
# 	direct = os.path.abspath(os.path.dirname(__file__))
# 	pro_dir = re.findall('(^.*)\/snpdb' , direct)[0]
# 	vcf_path = os.path.join(direct, 'vcf_files')
#
# 	#Collects all libraries to be compared.
# 	group_1_path = []
# 	for each in vcf1:
# 		vcf1_path = os.path.join(pro_dir, each)
# 		group_1_path.append(vcf1_path)
# 	group_2_path = []
# 	for each in vcf2:
# 		vcf2_path = os.path.join(pro_dir, each)
# 		group_2_path.append(vcf2_path)
# 	libraries = group_1_path + group_2_path
# 	count = len(libraries)
# 	libs = library_1 + library_2
# 	vcf_string = '_'.join(libs)
#
# 	#Determines the location of where analysis results will be stored.
# 	path = os.path.join(vcf_path, 'bcftools_isec_snpEff_%s_%s' % (vcf_string, datetime.datetime.utcnow().strftime("%Y-%m-%d")))
#
# 	#Checks to see if analysis has already been completed. If the analysis files are not present, bcftools is called.
# 	if os.path.isfile(path):
# 		pass
# 	else:
# 		#Checks to see if files have been zipped and indexed. Bcftools requires indexed vcf files.
# 		for fname in libraries:
# 			if os.path.isfile(fname):
# 				#zips and indexes vcf-files for bcftools
# 				try:
# 					subprocess.check_call(['bgzip', fname])
# 					subprocess.check_call(['tabix', '-p', 'vcf', '%s.gz' % fname])
# 				except IOError:
# 					pass
# 			elif os.path.isfile('%s.gz' % fname):
# 				print "files already zipped"
#
# 		#Creates a string of all zipped file
# 		zip_vcf = ''
# 		for each in libraries:
# 			zips = str(each) + '.gz'
# 			zip_vcf = zip_vcf + ' ' + zips
#
# 		#Runs the bcftools isec command to compare results. Files are outputed to three separate files.
# 		subprocess.check_call(["""bcftools isec -n -1 %s -p %s""" % (zip_vcf, path)],
# 		                      shell=True)
#
# 		#unzips files to return to the original state.
# 		for fname in libraries:
# 			try:
# 				subprocess.check_call(['gunzip', '%s.gz' % fname])
# 			except subprocess.CalledProcessError:
# 				print "File was not unzipped."
#
# 	#Opens the returned bcftools vcf files and counts the data.
# 	lib_effect = []
# 	total_2 = 0
# 	for i in range(0, count):
# 		vcf_reader = vcf.Reader(open ('%s/000%s.vcf' % (path, i), 'r'))
# 		lib = libs[i]
# 		high_effects = defaultdict(int)
# 		moderate_effects = defaultdict(int)
# 		modifier_effects = defaultdict(int)
# 		low_effects = defaultdict(int)
# 		total_counts = [0, 0, 0, 0, 0, '']
# 		for record in vcf_reader:
# 			effects = record.INFO['EFF']
# 			#Keeps track of what effect type each snp has. [high, moderate, low, modifier]
# 			impact_counts = [0, 0, 0, 0]
# 			#Places each type of impact into dictionary of the effect. SNPs with multiple impacts will have all impacts accounted for in the impact total.
# 			# i.e, SNPs with Downstream and Upstream effects will results in an addition to both impact counts.
# 			for x in effects:
# 				impact = x.split('(')[0]
# 				effect_list = x.split('(')[1]
# 				effect = effect_list.split('|')
# 				if effect[0] == "HIGH":
# 					impact_counts[0] += 1
# 					high_effects[impact] += 1
# 				elif effect[0] == "MODERATE":
# 					impact_counts[1] += 1
# 					moderate_effects[impact] += 1
# 				elif effect[0] == "MODIFIER":
# 					impact_counts[3] += 1
# 					modifier_effects[impact] += 1
# 				elif effect[0] == "LOW":
# 					impact_counts[2] += 1
# 					low_effects[impact] += 1
#
# 			#Counts the number of snps effected by each impact type. Snp is only counted once for each impact, i.e. if SNP has two modifying impacts, it is only counted once.
# 			total_counts[5] = lib
# 			total_counts[4] += 1
# 			if impact_counts[0] > 0:
# 				total_counts[0] += 1
# 			if impact_counts[1] > 0:
# 				total_counts[1] += 1
# 			if impact_counts[2] > 0:
# 				total_counts[2] += 1
# 			if impact_counts[3] > 0:
# 				total_counts[3] += 1
# 		lib_tuple = (dict(high_effects), dict(moderate_effects), dict(modifier_effects), dict(low_effects), total_counts)
# 		lib_effect.append(lib_tuple)
# 		total_2 += total_counts[4]
#
# 	return render_to_response('snpdb/test2.html', {"high_effects": dict(high_effects),
# 	                                               "moderate_effects": dict(moderate_effects),
# 	                                               "modifier_effects": dict(modifier_effects),
# 	                                               "low_effects": dict(low_effects),
# 	                                               "total_counts": total_counts,
# 	                                               "total_2": total_2,
# 	                                               "lib_effect": lib_effect,
# 	                                               "library1": library_1,
# 	                                               "library2": library_2,
# 	                                               "path": path}, context_instance=RequestContext(request))
#
