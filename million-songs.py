import random, time, math
import sys, os

range_left,range_right, output_F=sys.argv[1:]
range_left=int(range_left)
range_right=int(range_right)

def NUM_SONGS(file):
    songlist=dict()
    with open(file,"r") as f:
        for line in f:
            _,song,_=line.strip().split('\t')
            if song in songlist:
                songlist[song]+=1
            else:
                songlist[song]=1
    return songlist

def MAP_SONG_USERS(file,set_users=None, ratio=1.0):
    stu=dict()
    with open(file,"r") as f:
        for line in f:
            if random.random()<ratio:
                user,song,_=line.strip().split('\t')
                if not set_users or user in set_users:
                    if song in stu:
                        stu[song].add(user)
                    else:
                        stu[song]=set([user])
    return stu

def MAP_USER_SONG(file):
    list_user_songs=dict()
    with open(file,"r") as f:
        for line in f:
            user,song,_=line.strip().split('\t')
            if user in list_user_songs:
                list_user_songs[user].add(song)
            else:
                list_user_songs[user]=set([song])
    return list_user_songs

def Load_u(file):
    with open(file,"r") as f:
        u=map(lambda line: line.strip(),f.readlines())
    return u

def Index_s(file):
     with open(file,"r") as f:
         sti=dict(map(lambda line: line.strip().split(' '),f.readlines()))
     return sti

def Saving(r,songs_file,ofile):
    s2i=Index_s(songs_file)
    print "Saving ..."
    f=open(ofile,"w")
    for r_songs in r:
        indices=map(lambda s: s2i[s],r_songs)
        f.write(" ".join(indices)+"\n")
    f.close()

def theUnique(file):
    u=set()
    with open(file,"r") as f:
        for line in f:
            user,_,_=line.strip().split('\t')
            if user not in u:
                u.add(user)
    return u 

def songs_sorting(d):
    return sorted(d.keys(),key=lambda s:d[s],reverse=True)

def fl():
    sys.stdout.flush()

#l_rec: list of recommended songs
#u2s: mapping users to songs
#tau: 500
def AP(l_rec, sMu, tau):

    np=len(sMu)
    #print "np:", np
    nc=0.0
    mapr_user=0.0
    for j,s in enumerate(l_rec):
        if j>=tau:
            break
        if s in sMu:
        #print "s in sMu"
            nc+=1.0
            mapr_user+=nc/(j+1)
    mapr_user/=min(np,tau)
    return mapr_user

#l_users: list of users
#l_rec_songs: list of lists, recommended songs for users
#u2s: mapping users to songs
#tau: 500
def mAP(l_users, l_rec_songs, u2s, tau):
    mapr=0
    n_users=len(l_users)
    for i,l_rec in enumerate(l_rec_songs):
        if not l_users[i] in u2s:
            continue
        mapr+=AP(l_rec,u2s[l_users[i]], tau)
    return mapr/n_users

###
### PREDICTORS
###

class Pred:
    '''Implement generic predictor'''        
    
    def __init__(self):
        pass

    def Score(self,u2songs,  total_s):
        return {}

class Introduce_implement(Pred):
    '''Implement song-similarity based predictor'''

    def __init__(self, _intro_songs_users, EA=0, EQ=1):
        Pred.__init__(self)
        self.intro_songs_users = _intro_songs_users
        self.Q = EQ
        self.A = EA

   
    def Match(self,s,u_song):
        l1=len(self.intro_songs_users[s])
        l2=len(self.intro_songs_users[u_song])
        up = float(len(self.intro_songs_users[s]&self.intro_songs_users[u_song]))
        if up>0:
            dn = math.pow(l1,self.A)*math.pow(l2,(1.0-self.A))
            return up/dn
        return 0.0

    def Score(self,u2songs,  total_s):
        aggre_score={}
        for s in total_s:
            aggre_score[s]=0.0
            if not (s in self.intro_songs_users):
                continue
            for u_song in u2songs:
                if not (u_song in self.intro_songs_users):
                    continue
                s_match=self.Match(s,u_song)
                aggre_score[s]+=math.pow(s_match,self.Q)
        return aggre_score

class Introduce_implementc(Introduce_implement):
    '''Implement calibrated song-similarity based predictor''' 

    def __init__(self, _intro_songs_users, EA=0, EQ=1, f_hsongs=""):
        Introduce_implement.__init__(self, _intro_songs_users, EA, EQ)
        self.hsongs={}
        with open(f_hsongs,"r") as f:
            for line in f:
                s,v = line.strip().split()
                self.hsongs[s]=float(v)
        self.THETA = 0.5

    def select_theta(self,h):
        return self.THETA
        
    def calibrate(self, sco, song):
        h = self.hsongs[song]
        theta = self.select_theta(h)
        cali=sco
        if sco<h:
            cali = theta*sco/h
        elif sco>h:
            cali = theta+(1.0-theta)*(sco-h)/(1.0-h)
        return cali

    def Score(self, u2songs,  total_s):
        np = len(u2songs)
        aggre_score={}
        for s in total_s:
            aggre_score[s]=0.0
            if not (s in self.intro_songs_users):
                continue
            for u_song in u2songs:
                if not (u_song in self.intro_songs_users):
                    continue
                s_match=self.Match(s,u_song)
                aggre_score[s]+=math.pow(s_match,self.Q)/np
        for s in total_s:
            if s in self.hsongs:
                aggre_score[s]=self.calibrate(aggre_score[s],s)
            else:
                aggre_score[s]=0.0
        return aggre_score        
    
class User_similarity(Pred):

    '''Implement user-similarity based predictor'''
    
    def __init__(self, _u2s_tr, EA=0, EQ=1):
        Pred.__init__(self)
        self.u2s_tr = _u2s_tr
        self.Q = EQ
        self.A = EA
    
    

    def Score(self,u2songs,  total_s):
        aggre_score={}
        for u_tr in self.u2s_tr:
            if not u_tr in self.u2s_tr:
                continue
            w=float(len(self.u2s_tr[u_tr] & u2songs))
            if w > 0:
                l1=len(u2songs)
                l2=len(self.u2s_tr[u_tr])
                w/=(math.pow(l1,self.A)*(math.pow(l2,(1.0-self.A))))
                w=math.pow(w,self.Q)
            for s in self.u2s_tr[u_tr]:
                if s in aggre_score:
                    aggre_score[s]+=w
                else:
                    aggre_score[s]=w
        return aggre_score

###
### RECOMMENDERS
###

class Reco:

    '''Implements Recommender'''

    def __init__(self, _total_s):
        self.predictors=[]
        self.total_s=_total_s
        self.tau=500

    def Add(self,p):
        self.predictors.append(p)

class Liked_songs(Reco):

    '''Implements Aggregated Stochastic Recommender'''

    def __init__(self,_total_s):
        Reco.__init__(self,_total_s)
        self.Gamma=[]

    def Random_Index(self,n,distr):
        r=random.random()
        for i in range(n):
            if r<distr[i]:
                return i
            r-=distr[i]
        return 0
        
    def Random_recom(self,sorting, distr):
        nPreds=len(self.predictors)
        r=[]
        ii = [0]*nPreds
        while len(r)<self.tau:
            pi = self.Random_Index(nPreds,distr)
            s = sorting[pi][ii[pi]]
            if not s in r:
                r.append(s)
            ii[pi]+=1
        return r

    def songs_sorting(d):
        return sorted(d.keys(),key=lambda s:d[s],reverse=True)

    def Valid(self, T, users_te, calibrate, u2s_h, n_batch=10):
        aveEAP=0.0
        for t in range(T):
            rusers = users_te[t*n_batch:(t+1)*n_batch]
            rec=[]
            start=time.clock()
            for i,ru in enumerate(rusers):
                if ru in calibrate:
                    print "%d] scoring user %s with %d songs"%(i,ru,len(calibrate[ru]))
                else:
                    print "%d] scoring user %s with 0 songs"%(i,ru)
                fl()
                sorting=[]
                for p in self.predictors:
                    i_songs=[]
                    if ru in calibrate:
                        i_songs=songs_sorting(p.Score(calibrate[ru],self.total_s))
                    else:
                        i_songs=list(self.total_s)
                   
                    cleaned_songs = []
                    for x in i_songs:
                        if len(cleaned_songs)>=self.tau: 
                            break
                        if ru not in calibrate or x not in calibrate[ru]:
                             cleaned_songs.append(x)
                                            
                    sorting+= [cleaned_songs]
                    
                rec += [self.Random_recom(sorting, self.Gamma)]

            cti=time.clock()-start
            fl()
            # valuta la rec cn la map
            map_cur = mAP(rusers,rec,u2s_h,self.tau)
            aveEAP+=map_cur
            fl()
    
        print "Done!"

    def BasicReco(self, user, calibrate):
        sorting=[]
        for p in self.predictors:
            i_songs=[]
            if user in calibrate:
                i_songs=songs_sorting(p.Score(calibrate[user],self.total_s))
            else:
                i_songs=list(self.total_s)

            cleaned_songs = []
            for x in i_songs:
                if len(cleaned_songs)>=self.tau:
                    break
                if x not in calibrate[user]:
                    cleaned_songs.append(x)

            sorting += [cleaned_songs]

        return self.Random_recom(sorting, self.Gamma)

    def recommendation(self, l_users, calibrate):
        sti=time.clock()
        Liked_songs=[]
        for i,u in enumerate(l_users):
            if not (i+1)%10:
                if u in calibrate:
                    print "%d] %s w/ %d songs"%(i+1,l_users[i],len(calibrate[u])),
                else:
                    print "%d] %s w/ 0 songs"%(i+1,l_users[i]),
                fl()
            Liked_songs.append(self.BasicReco(u,calibrate))
            cti=time.clock()-sti
            if not (i+1)%10:
                print " tot secs: %f (%f)"%(cti,cti/(i+1))
        fl()
        return Liked_songs



print "range_left: %d , range_right: %d"%(range_left,range_right)
sys.stdout.flush()

train_file="train_triplets.txt"
evaluate_file="kaggle_visible_evaluation_triplets.txt"

print 'loading users in %s'%"kaggle_users.txt"
sys.stdout.flush()
users_v=list(Load_u("kaggle_users.txt"))

print 'default ordering by popularity'
sys.stdout.flush()
songs_ordered=songs_sorting(NUM_SONGS(train_file))

print  "loading unique users indexes"
uu=theUnique(train_file)
u2i = {}
for i,u in enumerate(uu):
    u2i[u]=i

print 'song to users on %s'%train_file
intro_songs_users=MAP_SONG_USERS(train_file)

print "converting users to indexes"
for s in intro_songs_users:
    s_filter = set()
    for u in intro_songs_users[s]:
        s_filter.add(u2i[u])
    intro_songs_users[s]=s_filter

del u2i

print 'user to songs on %s'%evaluate_file
calibrate=MAP_USER_SONG(evaluate_file)


EA = 0.15
EQ = 3

pr=Introduce_implement(intro_songs_users, EA, EQ)
instance = Liked_songs(songs_ordered)
instance.Add(pr)
instance.Gamma=[1.0]

r=instance.recommendation(users_v[range_left:range_right],calibrate)
Saving(r,"kaggle_songs.txt",output_F)