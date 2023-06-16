#############################################################################
# pyraflogin.py - Include your own routines in pyraf.
# To use this with pyraf it is BEST BEST BEST BEST to make an alias:
# alias pyraf="(export PYTHONSTARTUP=$HOME/iraf/pyraflogin.py; /home/albert/anaconda3/envs/iraf27/bin/pyraf)"
#
# This is a demo (and learning exercise for author) of ways to influence PyRAF/IRAF.
#
# uparm list format reminders...:
#
#   name,type,mode,default,min,max,prompt
#
#   type of parameter. Allowable values:
#      b : means parameter of boolean
#      i : means parameter of integer
#      r : means parameter of real
#      s : means parameter of string
#
#      f : for parameter of filename type. f may be followed by any combination of:
#         r  read access
#         w  write acces
#         e  presence of file (exists)
#         n absence of file.
#         Thus fw means test whether file given as a value of
#
#
# mode (allowable combinatins: a/h/q/hl/ql)
#   a auto (use the 'mode' value given.
#   h hidden no questions asked,
#   l learn  remember this when the task finishes
#   q query  always ask
#
# the default mode for each 'a' parameter.uu
# (insert (buffer-file-name))
#
# (ediff-current-file)
# (wg-python-fix-pdbrc)
# (find-file-other-frame "./.pdbrc")
# (wg-python-fix-pdbrc)   # PDB DASH DEBUG end-comments
#
# (setq mypdbcmd (concat (buffer-file-name) "<args...>"))
# (progn (wg-python-fix-pdbrc) (pdb mypdbcmd))
#
# (wg-astroconda-pdb)       # IRAF27
# (wg-astroconda3-pdb)      # CONDA Python3
#
# (set-background-color "light blue")
# (wg-python-toc)
#
# def fitserialize(files, start=1)
# 
# def newfits(filename,data,header=None):
#
# def newcube(filename,filelist,header=None):
#
# def shiftup(filename,count=1):
#
# def shiftside(filename,count=1):
#
# def our_figure():
#
# def our_surface(np2d,title="Surface Plot",xlabel="pixels",ylabel="pixels",zlabel="ADU",logflag=False):
#
# def our_histogram(mydata,title='Histogram'):
#
# def our_simple_plot(aline,x=None,title='Linear Plot',degree=None):
#
# def our_plot(aline,x=None,degree=None,legend=None):
#
# def our_plot_show(title='Linear Plot',usegrid=True):
#
# def our_scatter_plot(x,y,title=None,spkeywords={}):
#
# def our_plot3d(z):
#
# def pypar(task):                           #   irafIntrospect
#
# def pyparload(task,filename):                        #   irafIntrospect
#
# def pyfixbinning(imagename, newx,newy,xkwbin,ykwbin):
#
# def r2s (pobjra):
#
# def d2s (pobjdec):
#
# def s2r(rastr):  # PDB -DEBUG
#
# def s2d(decstr):
#
# def fixCCDSOFT(filename):
#
# def fakedispersion(filename,pos=466,width=10):
#
# def pyhelp(func):
#
# 2017-11-17T08:06:55-0700 wlg
# 2019-08-28T01:48:38-0600 wlg
#############################################################################
import os,sys
try:
   import numpy as np
   import re
   from pyraf      import iraf
   from astropy.io import fits
   import pprint
except:
   print >>sys.stderr,"Unable to find iraf. Running under anaconda?"
   sys.exit(1)


# permit some quick plots
try:
   from   mpl_toolkits.mplot3d import Axes3D
   import matplotlib.pyplot as plt
   from   matplotlib import cm
except exception:
   print("Could not find python matplotlib parts we wanted. Running under anaconda?")

pyraflogin_doc = """pyraflogin functions
fitserialize(files, start=1)
newfits(filename,data,header=None)
newcube(filename,filelist,header=None)
newmef(outname,filelist,modify=True)
shiftup(filename,count=1)
shiftside(filename,count=1)
our_figure()
our_surface(np2d,title="Surface Plot",
   xlabel="row",ylabel="column",zlabel="ADU",logflag=False)
our_histogram(mydata,title='Histogram',show=True,sigma=None)
our_simple_plot(aline,x=None,title='Linear Plot',degree=None)
our_plot(aline,x=None,degree=None,legend=None)
our_plot_show(title='Linear Plot',usegrid=True)
our_scatter_plot(x,y,title=None,spkeywords={})
our_plot3d(z)
pypar(task)
pyparload(task,filename)
pypath(task)
pyfixbinning(imagename, newx,newy,xkwbin,ykwbin)
r2s(pobjra)
d2s(pobjdec)
s2r(rastr)
s2d(decstr)
our_makemask(fitsname,limit=None,nstd=5)
fixCCDSOFT(filename)
fakedispersion(filename,pos=466,width=10)
fakeidentify(filename,pos=466,width=10)
pyhelp(func)
our_spectro_hint(sci,comp)
myhistory(filename)
"""


#############################################################################
# make sure this is a pyraf script. Don't want to do this
# for ordinary python.
# executable = sys.argv[0]
# while os.path.islink(executable):
#    executable = os.readlink(executable)
# if( os.path.split(executable)[1] != "pyraf"):
#    ... take some kind of action?
#
# Add some tasks. Not easy to deduce.
#    pypar an lpar that prints a dictionary for use as a **kwds for function.
#
#############################################################################


##############################################################################
# fitserialize Prepend a SerialNo. to the files.
#
##############################################################################
def fitserialize(files, start=1):
   """Prepend a serial number to the file."""
   dates = {}
   msgs = []
   try:
      for filename in files:
         with fits.open(filename) as f:
            h = f[0].header
            if('DATE-OBS' in h):
               dates.setdefault(filename,[]).append(h['DATE-OBS'])
   except Exception as e:
      print >> sys.stderr,"File %s failed, aborting." % filename + e.__str__()
      sys.exit(1)
   keys = list(dates.keys())
   keys.sort()
   for fname in keys:
      with fits.open(fname):
         h = f[0].header
         number = "1"+("%d" % start).zfill(4)
         h['serno'] = number
         start += 1
         os.rename(fname,number+'_'+fname)

# fitserialize

#############################################################################
# newfits - given a new "filename" and data with optional header
# make a new fits file. This handles hacking with numpy.arrays on the
# data of an existing file.
#############################################################################
def newfits(filename,data,header=None):
   """newfits(filename,data,header=None)
   Write a new fits file based on data. Will overwrite."""
   try:
      if(header is None):
         nf = fits.PrimaryHDU(data)
      else:
         nf = fits.PrimaryHDU(data,header)
      nf.writeto(filename,output_verify='fix',overwrite=True)
   except Exception as e:
      print >>sys.stderr,"newfits: error with operation %s" % e # soft fail.
      return
   print "# %s written !ds9 %s & " % (filename,filename)

# newfits

#############################################################################
# newcube - for the files in an array or at-file. open each file, use
# the first header; combine the files into a MEF
#############################################################################
def newcube(filename,filelist,header=None):
   """newcube(filename,filelist,header=None)
   Combine the data extents for files in 'filelist' (array or '@filename' string)
   and if the shapes are the same; and write, but do not clobber, the filename.
   Raise exceptions for bad conditions."""
   msg = ""
   if(type(filelist) == type('str') and filelist[0] == '@'):
      newlist = []
      with open(filelist[1:],'r') as flist:
         for fname in flist:
            if(not os.path.exists(fname)):
               newlist.append(fname.strip())
            else:
               msg += "Filelist %s -- filename %s not found.\n" % (filelist,fname)
      filelist = newlist
   elif (type(filelist) != type([])):
      msg += "Expecting array of file names or an @file."
   if(msg != ""):
      print >>sys.stderr,msg,"\nAborting"
      raise Exception(msg)
   mefdata = []
   shapes  = {}
   for fname in filelist:
      with fits.open(fname) as f:
         if(header == None):
            header = f[0].header           # impose param or first header as 'the' header
         d = f[0].data
         mefdata.append(d)
         shapekey = d.shape.__str__()
         if( shapekey not in shapes):
            shapes[shapekey] = []
         shapes[shapekey].append([fname])  # accumulate the fnames per shape.
   if(len(shapes) != 1):
      emsg = ""
      for k,v in shapes.items():
         emsg += "%s %s" % (k, '\n   '.join(v))
      raise Exception("Shapes differ: %s" % emsg)
   mef = np.array(mefdata)
   nf  = fits.PrimaryHDU(mefdata)
   nf.header = header
   nf.writeto(filename,output_verify='fix',overwrite=True)
   print "# %s written. !ds9 %s &" % (filename,filename)
   return mef                             # return the actual np array!
# newcube newcube('d30cube.fits','@l.l')

#############################################################################
# newmef -- merge several image/table files into one MEF.
#
#############################################################################
def newmef(outname,filelist,modify=True):
   """newmef(filelist) given a list of existing filenames, open each and
   add the data/header to a new mef file called outname.
newmef('testmef.fits',filelist)
ds9 -multiframe testmef.fits # open each HDU into its own frame
ds9 -mecube testmef.fits     # first HDU into one frame a cube slider for rest
ds9 testmef.fits             # open just first HDU in one frame
"""
   new_hdul = fits.HDUList()
   mainextname = None
   for i,fname in enumerate(filelist):
      with fits.open(fname) as f:
         for hdu in f:                   # e.g.: NICFPS files have 2 HDUs per file.
            if(modify):
               if('EXTNAME' in hdu.header): # guarantee a unique identifier.
                  if(mainextname is None):
                     extname = mainextname = hdu.header['EXTNAME']
                  else:
                     if(hdu.header['EXTNAME'] == mainextname):
                        extname = mainextname + "%d" % i
                  hdu.header['EXTNAME'] = extname
               else:
                  extname = 'EXTENSION'+ "%d" % i
                  hdu.header['EXTNAME'] = extname
            new_hdul.append(fits.ImageHDU(hdu.data,hdu.header))
   new_hdul.writeto(outname, overwrite=True)
# newmef newmef('testmef.fits',filelist)

#############################################################################
# shiftup
#
#############################################################################
def shiftup(filename,count=1):
   """shiftup(filename,count=1)
   Shift an image up by some count, then fill in the bottom
   with zeros. Used to expore image shifts."""
   f       = fits.open(filename)
   d       = f[0].data
   h       = f[0].header
   newdata = d
   newdata[count:,:] = d[0:-count,:]   # x/y reversed in python -> numpy.array
   newdata[0:count,:] = 0
   newname = 'x_' + filename
   newfits(newname,d,h)
   f.close()

# shiftup

#############################################################################
# shiftside
#
#############################################################################
def shiftside(filename,count=1):
   """shiftside(filename,count=1)
   Shift an image over by some count, then fill in the bottom
   with zeros. Used to expore image shifts.
   """
   with fits.open(filename) as f:
      f                  = fits.open(filename)
      d                  = f[0].data
      h                  = f[0].header
      newdata            = d
      newdata[:,count:]  = d[:,0:-count]   # x/y reversed in python -> numpy.array
      newdata[:,0:count] = 0
      newname = 'x_' + filename
      newfits(newname,d,h)

# shiftside

#############################################################################
# Hack plotting
# our_figure -- start a new figure
# Hack plotting starts with a call to our_figure() to re-initialize the
# plot functions. Then various calls to plotting calls, then finish with
# call to our_plot_show(). This is a consistent use of the plotting function.
# These hacks are for quick peeks. Use inline matplotlib calls for publication.
#
# our_surface
#############################################################################

_fig     = None  # manage an internal figure, reset with our_figure
_ax      = None

#############################################################################
# our_figure
#############################################################################
def our_figure():
   """our_figure()
   Reset our global fig and ax instances. Skip to allow overplotting."""
   global _fig,_ax,_legend
   _fig     = plt.figure(1,clear=True)
   _ax      = plt.axes()
   _legend  = []

# our_figure

##############################################################################
# our_surface
#
##############################################################################
def our_surface(np2d,title="Surface Plot",xlabel="row",ylabel="column",zlabel="ADU",logflag=False):
   """our_surface(np_2d_data,title="Surface Plot",xlabel="pixels",ylabel="pixels",
         zlabel="ADU",logflag=False)
   Our surface plotter"""
   our_figure()
   _ax      = _fig.gca(projection='3d')
   _ax.set_xlabel(xlabel)
   _ax.set_ylabel(ylabel)
   if( isinstance(np2d,str)):
      with fits.open(np2d) as f:
         img = f[0].data
   else:
      img = np2d
   if(logflag):
      _ax.set_zlabel(zlabel + ' [log10]')
      img[img <= 0] = 1.0
      z               = np.log10(img)
   else:
      _ax.set_zlabel(zlabel)
      z               = img
   xx,yy              = np.meshgrid(np.arange(img.shape[1]), np.arange(img.shape[0]))
   surf               = _ax.plot_surface(xx, yy, z, cmap=cm.Reds,
                               linewidth=0, antialiased=False)
   plt.title(title)
   plt.show()
# our_surface

def our_histogram(mydata,title='Histogram',show=True,sigma=None):
   """our_histogram(hist,title='Histogram',show=True,sigma=None)
   and hist = np.hist(mydata) # mydata.shape(x,)
   our_histogram(hist,title='Histogram',show=True) -- given a
   np.histogram(...) tuple and an optional title, and a boolean
   to go ahead andpop up a display. If sigma is offered, histogram
   is made for data > sigma of the mean.
   Uses our global _fig and _ax"""
   global _fig, _ax
   our_figure()
   pdata = mydata
   if(sigma is not None):
      sigma  = 7.0
      select = sigma*mydata.std()
      data   = mydata.flatten()
      pdata  = data[np.where(data > data.min() + sigma * data.std())]
      title  = "%s Selected > %6.2f σ (ADU > %7.3f)" % (title,sigma,select)
   #bins  = range(int(dmin),int(dmax),int((dmax-dmin)// 10))
   #myhist = np.histogram()
   #bins   = len(myhist[0])
   #print(f"myhist {myhist} {bins}")
   plt.title(title)
   _ax.hist(pdata)
   if(show):
      plt.show()

# our_histogram

#############################################################################
# our_plot - simple plotter for data
#
#############################################################################
def our_simple_plot(aline,x=None,title='Linear Plot',degree=None):
   """our_simple_plot(aline,x=None,title='Linear Plot')
   Reset
   our_figure(), and then Given 'aline' of a np.array of shape(n,), plot
   it. Use x if aline.shape == x.shape and is not None; if x
   is not provided or x=None -- then create and use
   tmpx=np.arange(x=1,aline.shape[0]+1).
   """
   ret = None
   global _ax,_fig
   our_figure()
   if(x is None):
      x = np.arange(1,aline.shape[0]+1)
   ignore = _ax.plot(x,aline,linewidth=0.5)
   if(degree is not None):
      minx,maxx    = x.min(),x.max()
      coefficients = np.polyfit(x,aline, int(degree))
      polynomial   = np.poly1d(coefficients)
      xs           = np.arange(minx,maxx, 0.1)
      ys           = polynomial(xs)
      ignore       = _ax.plot(xs,ys)
      ret = (ys[0],ys[-1])
   plt.title(title)
   plt.tight_layout()
   plt.show()
   return ret

# our_simple_plot

#############################################################################
# our_plot - simple plotter for data
#
#############################################################################
def our_plot(aline,x=None,degree=None,legend=None):
   """our_plot(aline,x=None,title='Linear Plot')
   Given 'aline' np.array of shape(n,), plot it. Use x if
   aline.shape == x.shape and is not None; if x=None then create arange
   x=1,aline.shape[0]+1. Use our _fig and _ax."""
   ret = None
   global _ax,_fig,_legend
   if(x is None):
      x = np.arange(1,aline.shape[0]+1)
   _legend.append(legend)
   ignore = _ax.plot(x,aline,linewidth=0.5,label=legend)
   if(degree is not None):
      minx,maxx    = x.min(),x.max()
      coefficients = np.polyfit(x,aline, int(degree))
      polynomial   = np.poly1d(coefficients)
      xs           = np.arange(minx,maxx, 0.1)
      ys           = polynomial(xs)
      ignore = _ax.plot(xs,ys,label=legend)
      ret = (ys[0],ys[-1])
   return ret

# our_plot

#############################################################################
# our_plot_show - Show the plot with a title
#############################################################################
def our_plot_show(title='Linear Plot',usegrid=True):
   """our_plot_show(title='Linear Plot',usegrid=True
   After our_* show the plot."""
   global _fig,_ax,_legend
   leg = [ [l,""][l is None] for l in _legend]
   _ax.legend(leg)
   if(usegrid):
      _ax.grid()
   plt.title(title)
   plt.tight_layout()
   print "# %s " % leg
   plt.show()

# our_plot_show


#############################################################################
# our_scatter_plot
#############################################################################
def our_scatter_plot(x,y,title=None,spkeywords={}):
   """our_scatter_plot(x,y,title=none)"""
   global _ax,_fig
   ignore = _ax.scatter(x,y,**spkeywords)
   plt.title(title)
   plt.tight_layout()
# our_scatter_plot


#############################################################################
# our_polt3d
#############################################################################
def our_plot3d(z):
   """our_plot3d(z)
     Plot an image: our_plot3d(z) where z.shape (x,y) 2d array """
   fig   = plt.figure()
   ax    = fig.gca( projection='3d')
   X,Y   = np.mgrid[0:z.shape[0],0:z.shape[1]]
   ignore  = ax.plot_surface(X,Y, z, rstride=1, cstride=1,cmap=plt.cm.gray, linewidth=0)
# our_plot3d(d) our_plot3d(d[1:100,1:100])


#############################################################################
# Steps:
# 1) use iraf.osfn and makesure a param file is available
# 2) write the task
# 3) tie to iraf with iraf.IrafTaskFactory

#############################################################################
# guarantee pypar.par on first use. Use single ticks for quotation.
# Force a default uparm file if not found.
#############################################################################
pyparfname = iraf.osfn("uparm$pypar.par")  # need the name later

if(not os.path.islink(pyparfname)):        # force the parm file if first time
   with open(pyparfname,'w') as f:
      print >>f,"""funcname,s,a,'imstat',,,'iraf function for introspection'
mode,s,h,'al'
"""

#############################################################################
# pypar - the actual text of the function into the namespace
#############################################################################
def pypar(task):                           #   irafIntrospect
   """Print a dictionary form of lpars"""
   m     = globals()['iraf']
   ifunc = getattr(m,task)
   pp    = pprint.PrettyPrinter(indent=4)
   d     = {}
   for p in ifunc.getParList():
      v = p.__dict__['value']
      d[p.__dict__['name']] = v
   pp.pprint(d)  # print d
# pypar

# tie pypar into the mix as a "task"
iraf.IrafTaskFactory(prefix='',taskname='pypar', suffix='', pkgbinary=None,
       value=pyparfname, function=pypar)

#############################################################################
################################## pyparload ################################
#############################################################################

#############################################################################
# guarantee pyparload.par on first use. Use single ticks for quotation.
# Force a default uparm file if not found.
#############################################################################
pyparloadfname = iraf.osfn("uparm$pyparload.par")  # need the name later

if(not os.path.islink(pyparloadfname)):        # force the parm file if first time
   with open(pyparloadfname,'w') as f:
      print >>f,"""funcname,s,a,'imstat',,,'iraf function for introspection'
mode,s,h,'al'
"""

#############################################################################
# pyparload - the actual text of the function into the namespace
#############################################################################
def pyparload(task,filename):                        #   irafIntrospect
   """From pyraf documentation. Print a dictionary form of lpars"""

   m     = globals()['iraf']
   try:                                              # get the task or die
      ifunc = getattr(m,task)
   except Exception as e:
      print >>sys.stderr,"task %s" % task, " not found"
      return                                         # fail  quietly

   try:
      with open(filename,'r') as f:                  # get the potential pypar dict.
         txt = f.read()
   except Exception as e:
      print >>sys.stderr,"Issue opening filename %s" % filename
      return

   try:
      newdict = eval(txt)                            # make new dict or die
   except Exception as e:
      print >>sys.stderr,"pyparload: %s" % filename, " error %s" % e
      return                                         # fail  quietly

   msg   = ""                                        # initial error message
   comma = ""                                        # comma separated list
   for k,v in newdict.items():
      if(len(ifunc.getAllMatches(k)) == 0):
         msg += "%s%s" % (comma,k)                   # oops a user's key not in
         comma = ", "                                # task's list add errors
   if( msg != ""):
      print >>sys.stderr,"parameters %s not found" % msg
      return                                         # fail  quietly

   for p in ifunc.getParList():                      # now scan the tasks
      n = p.__dict__['name']                         # see if task name missing
      if(n not in newdict):
         msg += "%s%s" % (comma,n)                   # add to errors

   if(msg != ""):                                    # if any msgs issue them
      print >>sys.stderr,"Task %s does not have parameters: %s" % (task,msg)
      return                                         # fail quietly

   # names supplied match names wanted
   for p in ifunc.getParList():
      n = p.__dict__['name']
      v = newdict[n]
      p.__dict__['value'] = v                       # duck-type add

# pyparload

def pypath(task):
   """Show the python path"""
   for p in sys.path:
      print(p)
# pypath

# tie pyparload into the mix as a "task"
iraf.IrafTaskFactory(prefix='',taskname='pyparload', suffix='', pkgbinary=None,
   value=pyparloadfname, function=pyparload)

iraf.IrafTaskFactory(prefix='',taskname='pypath', suffix='', pkgbinary=None,
   value="", function=pypath)


#############################################################################
############################## DONE WITH PYPAR ##############################
#############################################################################

#############################################################################
# Guarantee parameter file for pyfixbinning task, force a default uparm file
# if not found.
#############################################################################
binname    = iraf.osfn("uparm$pyfixbinning.par")
if(not os.path.islink(binname)):
   with open(binname,'w') as f:
      print >>f,"""imagename,s,a,'default.fits',,,'image name to rebin'"""
      print >>f,"""newx,i,a,2,1,5,'New XBINNING factor'"""
      print >>f,"""newy,i,a,2,1,5,'New YBINNING factor'"""
      print >>f,"""xkwbin,s,a,'XBINNING',,,'Keyword for x binning'"""
      print >>f,"""ykwbin,s,a,'YBINNING',,,'Keyword for y binning'"""
      print >>f,"""mode,s,h,'ql'"""
# inline

#############################################################################
# This task is not quite fully formed. It only takes one file at a time
# and updates in place.
# TODO: Add the @innames and x_//@outname trick
#       Multiple comma separated names etc.
#############################################################################
def pyfixbinning(imagename, newx,newy,xkwbin,ykwbin):
   """pyfixbinning(imagename, newx,newy,xkwbin,ykwbin)
   The parameters in order by name in the uparm file.
   Only 2x2 and 3x3 are supported. The binning has to be symmetric.
   """
   try:
      msg      = "imagename"                         # breadcrumb for each thing that can blow up.
      filename = imagename
      msg      = "newx"
      pnewx    = int(newx)                           # make sure thisis an integer
      msg      = "newy"
      pnewy    = int(newy)                           # make sure thisis an integer
      msg      = "xkwbin"
      xkeyword = xkwbin
      msg      = "ykwbin"
      ykeyword = ykwbin
   except Exception as e:
      # Oops cant do our thing.
      print >>sys.stderr,"pyfixbinning: Bad parameter(s) %s\n%s" % (msg,e)
      raise

   if(pnewx != pnwey and pnewx not in [2,3]):
      raise Exception("X spanning for 2x2 or 3x3 supported, You asked for %dx%d" % (pnewx,pnewy))

   with fits.open(filename, mode='update') as f:
      d = f[0].data
      h = f[0].header

      if(type(xkwbin) != type("") or type(ykwbin) != type("")): # dont duck-type here.
         raise Exception("Binning keywords are not strings.")
      if(xkwbin not in h or ykwbin not in h):
         raise Exception("Binning keywords not found, you specified X=%s, Y=%s" % (xkwbin,ykwbin))

      # hack other elif's for 4,5, etc. if you need them.
      if(pnewx == 3):
         dd=d[0:-2:3, 0:-2:3]   + d[1:-1:3,1:-1:3]  + d[2::3,2::3]
      elif(pnewx == 2):
         dd=d[0:-2:2, 0:-2:2]     + d[2::2,2::2]

      f[0].data            = dd                         # put things back, ready for update
      n1,n2                = dd.shape                   # get the new shape
#      f[0].header.set['NAXIS1']   = n2
#      f[0].header.set['NAXIS2']   = n1
#      f[0].header.set('XBINNING') = pnewx
#      f[0].header.set('YBINNING') = pnewy
      # the default close will update in place the file.

# pyfixbinning

#############################################################################
# tie pyfixbinning into the mix as a "task"
#############################################################################
iraf.IrafTaskFactory(prefix='',taskname='pyfixbinning', suffix='', pkgbinary=None,
    value=binname, function=pyfixbinning)


#############################################################################
########################## DONE WITH PYFIXBINNING ###########################
#############################################################################

##############################################################################
# r2s Convert RA in degrees to sexa in hours
#
##############################################################################
def r2s(pobjra):
   """r2s(pobjra)
   Convert RA in degrees (float or string) to sexigesimal in hours.
   Leading zero for pretty pringing. This is related to a custom
   postgresql function. Negative angles for hour angles can be supported
   by this code.
   Should test for -24 < pobjra < 24
   fmt:  '-HH:MM:SS.ss' ' HH:MM:SS.ss' ftmlen = %12s
   """

   pobjra = float(pobjra) / 15.0;         #  divide degrees into hours first

   idegs = math.floor(pobjra);            #  get the degrees part
   isecs = 60.0 * (pobjra - idegs);       #  get the
   imins = math.floor(isecs);             #  Nail down minutes
   isecs = 60.0 * (isecs - imins);        #  get the float seconds
   sign  = " "

   if(idegs < 0):
      pdegs = "%d" % (0.0 - idegs)
      sign  = "-"
   elif (idegs >= 0 and idegs < 10):
      pdegs = "0%d" % idegs
   else:
      pdegs = "%d" % idegs
   pdegs = sign + pdegs

   pmins = '%d' % imins
   if imins < 10:
      pmins = '0' + pmins

   psecs = "%.2f" % isecs
   if isecs < 10.0:
      psecs = '0' + psecs

   ret = "%s:%s:%s" % (pdegs, pmins, psecs)

   return ret;

# r2s

##############################################################################
# d2s - convert declination in degrees to sexadecimap
#
##############################################################################
def d2s(pobjdec):
   """d2s(pobjdec)
   Convert declinatin in degrees (float or string) to sexigecimal in
   degrees.  Leading zero for pretty pringing. This is related to a
   custom postgresql function.

   """
   psign = '+';
   pobjdec = float(pobjdec)
   if  pobjdec < 0.0:
      psign='-' ;
      pobjdec = 0.0 - pobjdec;

   idegs = math.floor(pobjdec);            #  get the degrees part
   isecs = 60.0 * (pobjdec - idegs);  #  borrow psec, to get the minutes
   imins = math.floor(isecs);              #  Nail down minutes
   isecs = 60.0 * (isecs - imins);    #  get the float seconds

   pdegs = "%d" % idegs
   if idegs < 10.0:
      pdegs = '0' + pdegs        #  pad degrees with leading zero

   pmins = "%d" % imins
   if imins < 10:
      pmins = '0' + pmins;        #  pad minutes with leading zero

   psecs = "%.2f" % isecs
   if isecs < 10.0:
      psecs = '0' +  psecs       #  pad seconds with leading zero

   ret = "%s%s:%s:%s" % (psign, pdegs, pmins, psecs)

   return ret;

# d2s

##############################################################################
# s2r - convert a sexigesimal RA TO a floating point degrees.
#   input is string hh:mm[:ss.s] Will take ra.ddddd ra:mm.mmm as well.
#   The truncated forms appear in SIMBAD query output.
##############################################################################
coordre = re.compile(r'[dmsh: ]+')
def s2r(rastr):  # PDB -DEBUG
   """s2r(rastr)
   Convert a sexigesimal RA to a floating point degrees."""
   if(type(rastr) != type(1.0)):  # PDB -DEBUG
      ra = None
      try:
         parts = list(map(float, (coordre.split(rastr) + 3 * ['0']) [:3]))
         if(len(parts) == 1):
            ret = parts[0] * 15.0     # a straight hrs.xxxxxxx
         elif(len(parts) == 2):
            ret = (parts[0] + (parts[1] / 60.0)) * 15.0
         elif(len(parts) >= 3):
            parts = parts[:3]   # guarantee guarantee 3 parts.
            ret = (parts[0] + (parts[1] / 60.0) + parts[2]/3600.0) * 15.0
      except:
         raise ValueError("s2r: Unable to convert %s to degrees ra" % rastr)
   else:
      ret = rastr
   return ret

# def s2r  s2r('12.34567') s2r('12:13.4567') s2r('12:13:45.67') s2r('20 08 24')

##############################################################################
# s2d - convert a sexigesimal Dec TO a floating point degrees.
#   input is string hh:mm[:ss.s] Will take ra.ddddd ra:mm.mmm as well.
#   The truncated forms appear in SIMBAD query output.
##############################################################################
def s2d(decstr):
   """Convert a sexigesimal Dec to a floating point degrees. IRAF upon
   occasion will return a declination in degrees (float or string)"""
   if(type(decstr) != type(1.0)):
      sign = 1.0
      try:
         parts = list(map(float, coordre.split(decstr)))
         if(parts[0] < 0.0):
            sign = -1.0
            parts[0] = -1 * parts[0]
         if(len(parts) == 1):
            dec = parts[0]     # a straight d.xxxxxxx
         elif(len(parts) == 2): # d m.xxxxxx
            dec = (parts[0] + (parts[1] / 60.0))
         elif(len(parts) == 3): # d m s.xxxx
            dec = (parts[0] + (parts[1] / 60.0) + parts[2]/3600.0)
         ret = sign*dec
      except:
         raise ValueError("r2d: Unable to convert %s to degrees dec" % decstr)
   else:
      ret = decstr

   return ret

# s2d s2d('-12.34567') s2d('12:34:56.7') s2d('-12:34:56.7') s2d('-12 34 56.7')

#############################################################################
# make a bad pixel mask called mask.pl and a mask.fits file, where good
# pixels are zero and bad pixels their original value in fitsname image.
#############################################################################
def our_makemask(fitsname,limit=None,nstd=5):
   """Make a mask.pl file where px > val in FITS fitsname.
      Also create a mask.fits file where good pixels are zero and
      bad pixels are their original value in fitsname image.
      Optional:  limit is a base mean. If None, then
      limit is computed as the mean + (nstd * stddev) for the image.
   """
   with fits.open(fitsname) as f:
      d      = f[0].data
   if(limit is None):
      mean   = d.mean()         # 4.3325748
      stddev = d.std()          # 249.91281
      limit  = mean + (nstd * stddev)
   bads = np.where(d > limit)
   if(len(bads[0]) > 0):
      with open('mask.pl','w') as o:
         for c,l in zip(*bads):
            print >>o,l+1,c+1  # 1's based and in FORTRAN index order
      img       = np.zeros_like(d)
      img[bads] = d[bads]
      newfits('mask.fits',img)
      msg = "# Parameters: %d %7.4f %7.4f %7.4f " % (len(bads[0]),mean,stddev,limit)
   else:
      msg = "No bad pixels at limit %f found." % (limit)
   return bads,
# our_makemask

#############################################################################
# Fix CCDSoft keywords fubar, fake a dispersion along a cal image for tracing.
#############################################################################
def fixCCDSOFT(filename):
   """fixCCDSOFT(filename)
   Fix CCDSOFT bad FITS headers"""
   f = fits.open(filename,mode='update')
   h    = f[0].header
   gain = h['E-GAIN']  # replace the E-GAIN with proper GAIN KW.
   del h['E-GAIN']
   h['GAIN'] = gain    # change the keyword.
   f.close(output_verify='silentfix') # blank pad, case for exponential.
# fixCCDSOFT

##############################################################################
# fakedispersion
#
##############################################################################
_fakedisparray = np.array([0.0, 0.0, 0.0, 0.0149098, 0.20077242, 0.61701065, 1.0,
                           0.61701065, 0.20077242, 0.0149098, 0.0, 0.0, 0.0])
def fakedispersion(filename,pos=466,width=10):
   """fakedispersion(filename,pos=466,width=13)
   Given <filename> put a fake Gaussian dispersion across column pos,
   with width and write sci_<filename> as output. Will overwrite filename"""
   try:
      if(width > 13):   # max set by _fakedisparray array of normalized gaussian vals.
         width = 13
      filename = filename.strip()
      msg             = "Open file"
      f               = fits.open(filename)
      msg             += '..get data'
      h               = f[0].header
      d               = f[0].data
      msg             += '..do work'
      __fakelen       = _fakedisparray.shape[0]
      __fakedisparray = _fakedisparray[__fakelen//2 - width//2: 1+__fakelen//2 + width//2]
      width           = __fakedisparray.shape[0]
      __fakedisparray = __fakedisparray.reshape(width,1)
      if(iraf.onedspec.dispaxis == 1):
         v = d[pos-width//2:1+pos+width//2,:]
         maxscale = v.max()//2
         d[pos-width//2:1+pos+width//2,:] = d[pos-width//2:1+pos+width//2,:] + (
            d[pos-width//2:1+pos+width//2,:] + (__fakedisparray * maxscale))
      else:
         v = d[:,pos-width//2:1+pos+width//2]
         maxscale = v.max()//2
         d[:,pos-width//2:1+pos+width//2] = d[:,pos-width//2:1+pos+width//2] + (
            d[:,pos-width//2:1+pos+width//2] ( __fakedisparray * maxscale))
      newname = 'sci_'+filename
      print "# Writing %s" % newname
      nf = fits.PrimaryHDU(d)
      msg += '..set header'
      nf.header = h
      msg += '..writeto'
      nf.writeto(newname,output_verify='fix',overwrite=True)
   except Exception as e:
      print  >>sys.stderr,"fakedispersion('{:s}',{:d},{:d})".format(filename, pos, width)
      print  >>sys.stderr,"fakedispersion error: file: '{:s}' msg: {:s}".format(filename, msg)
      print  >>sys.stderr,"position {:d}, width {:d}, shape {}".format(pos,width,d.shape)
      print  >>sys.stderr,"NAXIS1 = {:d}, NAXIS2 {:d}".format(h['NAXIS1'],h['NAXIS2'])
      print  >>sys.stderr,e.__str__()
# fakedispersion fakedispersion('190920_sun_1.fits',1075,7)

def fakeidentify(filename,pos=466,width=10):
   """fakeidentify(filename,pos=466,width=7)
   Given <filename> put a fake identify dispersion across column pos,
   with width and write sci_<filename> as output. Will overwrite filename"""
   try:
      if(width > 13):
         width = 13
      msg             = "Open file"
      f               = fits.open(filename)
      msg             += '..get data'
      h               = f[0].header
      d               = f[0].data
      msg             += '..do work'
      __fakelen       = _fakedisparray.shape[0]
      __fakedisparray = _fakedisparray[__fakelen//2 - width//2: 1+__fakelen//2 + width//2]
      width           = __fakedisparray.shape[0]
      __fakedisparray = __fakedisparray.reshape(width,1)
      if(iraf.onedspec.dispaxis == 1):
         v = d[pos-width//2:1+pos+width//2,:]
         d[pos-width//2:1+pos+width//2,:] = d[pos-width//2:1+pos+width//2,:] + (
            d[pos-width//2:1+pos+width//2,:] * __fakedisparray )
      else:
         v = d[:,pos-width//2:1+pos+width//2]
         d[:,pos-width//2:1+pos+width//2] = d[:,pos-width//2:1+pos+width//2] + (
            d[:,pos-width//2:1+pos+width//2] * __fakedisparray )
      newname = 'sci_'+filename
      print "# Writing %s" % newname
      nf = fits.PrimaryHDU(d)
      msg += '..set header'
      nf.header = h
      msg += '..writeto'
      nf.writeto(newname,output_verify='fix',overwrite=True)
   except Exception as e:
      print >>sys.stderr,"fakeidentify error. {:s} {:s}".format(filename, msg)
      print >>sys.stderr,"position {:d}, width {:d}, shape {}".format(pos,width,d.shape)
      print >>sys.stderr,"NAXIS1 = {:d}, NAXIS2 {:d}".format(h['NAXIS1'],h['NAXIS2'])
      print >>sys.stderr,e
# fakeidentify fakeidentify('190920_sun_1.fits',1075,7)


def pyhelp(func=None):
   """Given a fundtion, if func.__doc__ exists, print it."""
   if(func is None):
      print pyraflogin_doc
   else:
      try:
         print func.__doc__
      except:
         print "No help available"
# fakedispersion

def our_spectro_hint(sci,comp):
   """Hash up hints for reduction"""
   sciname = '.'.join(sci.split('.')[:-1])
   compname = '.'.join(comp.split('.')[:-1])
   s1 = """apall input=%s output=%s format="multispec" b_niterate = 7  b_sample = "-15:-7,7:15" t_function = "chebyshev" t_order = 7 t_niterate = 3 clean+""" % (sci, sci)
   print s1


# our_spectro_hint

def myhistory(filename):
   """Copy the history into provided filename"""
   import readline
   hist = [ readline.get_history_item(i)
             for i in xrange(1, readline.get_current_history_length() + 1)
             ]
   if(os.path.exists(filename)):
      with open(filename,'w') as f:
         f.write('\n'.join(hist))
      print("# Wrote history to {}".format(filename))
# myhistory

#############################################################################
#                      REPORT FILE HAS BEEN LOADED
#############################################################################
lfn = iraf.osfn("home$/pyraflogin.py")
print "\n%s loaded.\n " % lfn
print "home/pyraflogin: newfits     pypar    pyfixbinning   pyparload shiftup "
print "                 r2s         d2s      s2r            s2d"
print "                 fixCCDSOFT  fakedispersion"
print "                 our_surface our_histogram our_plot our_plot3d"
print " pyhelp(func) - safely print func.__doc__ if exists."
print "pyraflogin loaded"
