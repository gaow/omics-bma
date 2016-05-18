GSL_VERSION := 1.16
MOSEK_VERSION := 7.1.0.44
swig_opts := -c++ -python -O -shadow -keyword -w-511 -w-509

install:
	python setup.py install

install_libs: igsl ieqtlbma ideepdish

igsl:
	cd external/gsl && \
	./configure CFLAGS="-O3 -fPIC" --prefix=`pwd` && make && make install

ieqtlbma:
	swig $(swig_opts) -o src/pyeqtlbma/pyeqtlbma_wrap.cxx src/pyeqtlbma/pyeqtlbma.i && \
	mv src/pyeqtlbma/pyeqtlbma.py src/pyeqtlbma/__init__.py

ideepdish:
	cd external/deepdish && \
	python setup.py install

download: eqtlbma gsl deepdish mosek

eqtlbma:
	cd external && \
	git clone -b omics-bma https://github.com/gaow/eqtlbma.git

gsl:
	cd external && \
	rm -rf gsl* && \
	wget ftp://ftp.gnu.org/gnu/gsl/gsl-$(GSL_VERSION).tar.gz && \
	tar -zxvf gsl-$(GSL_VERSION).tar.gz && \
	ln -s gsl-$(GSL_VERSION) gsl

mosek:
	cd external && \
	rm -f mosektoolslinux64x86.tar.bz2 && \
	wget http://download.mosek.com/stable/$(MOSEK_VERSION)/mosektoolslinux64x86.tar.bz2 && \
	tar jxvf mosektoolslinux64x86.tar.bz2

deepdish:
	cd external && \
	git clone https://github.com/gaow/deepdish.git

clean:
	rm -rf build dist *.egg-info

.PHONY: eqtlbma gsl deepdish mosek install_libs igsl ieqtlbma ideepdish clean
