import re
import XRootD.client

def ls(server, directory):
  status, dirlist = XRootD.client.FileSystem(server).dirlist(directory)
  if not status.ok:
    raise Exception(
        'Got %d when requesting the contents of root://%s%s: %s' \
        % (status.code, server, directory, status.message)
    )
  return dirlist

def glob(server, directory, regex):
  return [
      'root://' + server + directory + '/' + f.name
      for f in ls(server, directory) if re.search(regex, f.name) 
  ]

def dirname(sample, datetime):
  return '/data/phedex/user/ezhemchu/crab/%(sample)s/%(sample)s/crab_%(sample)s/%(datetime)s' \
         % { 'sample': sample, 'datetime': datetime }

def file_job_id(filename):
  return int(re.search(r'(\d+)\.root$', filename).group(1))

def files_re(sample, datetime):
  server = 'cms-phedex.lxfarm.mephi.ru'
  directory = dirname(sample, datetime)
  files = []
  for d in ls(server, directory):
    files += glob(server, '%s/%s' % (directory, d.name), r'GEN-SIM-RAW-RECO_\d+\.root$')
  files.sort(key = file_job_id)
  return files

def files_range(sample, datetime, max):
  directory = dirname(sample, datetime)
  return [
      'root://cms-phedex.lxfarm.mephi.ru%s/%04d/%s_GEN-SIM-RAW-RECO_%d.root' \
      % (directory, i / 1000, sample, i) \
      for i in range(1, max + 1)
  ]

# assumes sets are sorted by file_job_id
def merge(sets):
  result = []
  i = [ 0 for s in sets ]
  n = len(sets)
  jid = 0
  while True:
    s = 0
    while s < n and i[s] >= len(sets[s]): s += 1
    if s == n: break
    min_s = s
    min_f = sets[s][i[s]]
    min_j = file_job_id(min_f)
    while s < n:
      if i[s] < len(sets[s]):
        j = file_job_id(sets[s][i[s]])
        if j < min_j:
          min_s = s
          min_f = sets[s][i[s]]
          min_j = j
      s += 1
    if min_j > jid:
      result.append(min_f)
      jid = min_j
    i[min_s] += 1
  return result

# W -> tau nu
# 37k events
def WToTauNu_190830():
  return files_re('WToTauNu_13TeV', '190830_120718')

# Z -> tau tau
# 50k events
def ZToTauTau_191020():
  return files_range('ZToTauTau_13TeV', '191020_113630', 50)

# W_R -> tau_R nu_R
# 43k events
def WRToTauNu_191023():
  return merge([
      files_re('WRToTauNu_13TeV', '191021_135640'),
      files_re('WRToTauNu_13TeV', '191021_140726'),
      files_re('WRToTauNu_13TeV', '191023_114855')
  ])

### cmsRun complains about duplicate event numbers. If we change Run or Lumi,
### we can have extra 31k events
## W_R -> tau_R nu_R
## 75k events
## (accidentally submitted twice in the first run)
#def WRToTauNu_191023():
#  return files_re('WRToTauNu_13TeV', '191021_135640') \
#       + files_re('WRToTauNu_13TeV', '191021_140726') \
#       + files_re('WRToTauNu_13TeV', '191023_114855')

# W -> tau nu
# tau decays only into hadrons
# tau is made unpolarized with Tauola
# 47k events
def WToTauNu_hadr_nopolar_191031():
  return files_re('WToTauNu_hadr_nopolar_Tauola', '191031_182731')

# W -> tau nu
# tau decays only into hadrons
# tau is made right-handed with Pythia
# 29k events
def WToTauNu_hadr_right_191102():
  # the files were erroneously named 'nopolar'
  return files_re('WToTauNu_hadr_nopolar_Pythia', '191102_114058')

# W -> tau nu
# tau decays only to hadrons
# tau is made unpolarized with Tauola
# tau pt is forced to be larger than 50 GeV by Pythia
# 50k events
def WToTauNu_hadr_nopolar_191117():
  return files_range('WToTauNu_hadr_nopolar_Tauola', '191117_124910', 50)

# W -> tau nu
# tau decays only to hadrons
# tau is left-handed (SM)
# tau pt is forced to be larger than 50 GeV by Pythia
# 921k events
def WToTauNu_191120():
  return files_re('WToTauNu', '191120_122523')

# W -> tau nu
# tau decays only to hadrons
# tau is right-handed (Pythia)
# 965k events
def WToTauNu_hadr_right_191127():
  # the files were erroneously named 'nopolar'
  return files_re('WToTauNu_hadr_nopolar_Pythia', '191127_113834')

# W -> ta u nu
# tau decays only to hadrons
# tau is right-handed (Pythia)
# tau pt > 50 GeV (Pythia)
# 904k events
def WToTauNu_hadr_right_191202():
  return files_re('WToTauNu_hadr_right_Pythia', '191202_124015')
