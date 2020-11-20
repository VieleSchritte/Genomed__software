from django.conf.urls import url

from . import views

app_name = 'cognation'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^calculate/$', views.calculate, name='calculate'),
    url(r'^save_allele/', views.save_allele, name='save_allele')
]
