from django.conf.urls.defaults import *
from views import hello, current_datetime, hours_ahead
from pwms import views
import settings

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    ('^$', views.search),
    ('^hello/$', hello),
    ('^time/$', current_datetime),
    (r'^time/plus/(\d{1,2})/$', hours_ahead),
    (r'^search-form/$', views.search_form),
    (r'^search/$', views.search),
    (r'^download/$', views.download),
    (r'^about/$', views.about),
    (r'^references/$', views.references),
    (r'^browse/$', views.browse),
    (r'^matrix/$', views.matrix),
    (r'^details/Y(.*)/$', views.systematic),
    (r'^details/(.*)/$', views.details),
    (r'^CHIP/(.*)/$', views.genes),
    (r'^KO/(.*)/$', views.genes),
    (r'^Pairs/(.*)/$', views.pairs),
    (r'^example/(.*)/$', views.example),
    (r'^static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.STATIC_ROOT}),
    (r'^templates/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.TEMPLATE_ROOT}),
    # Example:
    # (r'^TFDB/', include('TFDB.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # (r'^admin/', include(admin.site.urls)),
)
