import re
import XRootD.client

def glob(server, directory, regex):
  status, dirlist = XRootD.client.FileSystem(server).dirlist(directory)
  if not status.ok:
    raise Exception(
        'Got %d when requesting the contents of root://%s%s: %s' \
        % (status.code, server, directory, status.message)
    )

  return [
      'root://' + server + directory + '/' + f.name
      for f in dirlist if re.search(regex, f.name) 
  ]

def dirname(sample, datetime):
  return '/data/phedex/user/ezhemchu/crab/%(sample)s/%(sample)s/crab_%(sample)s/%(datetime)s/0000' \
         % { 'sample': sample, 'datetime': datetime }

def files_re(sample, datetime):
  files = glob(
      'cms-phedex.lxfarm.mephi.ru',
      dirname(sample, datetime),
      r'GEN-SIM-RAW-RECO_\d+\.root$'
  )
  files.sort(key = lambda f: int(re.search(r'(\d+)\.root$', f).group(1)))
  return files

def files_range(sample, datetime, max):
  directory = dirname(sample, datetime)
  return [
      'root://cms-phedex.lxfarm.mephi.ru%s/%s_GEN-SIM-RAW-RECO_%d.root' \
      % (directory, sample, i) \
      for i in range(1, max + 1)
  ]

# W -> tau nu
# 37k events
def WToTauNu_190830():
  return files_re('WToTauNu_13TeV', '190830_120718')

# Z -> tau tau
# 50k events
def ZToTauTau_191020():
  return files_range('ZToTauTau_13TeV', '191020_113630', 50)
