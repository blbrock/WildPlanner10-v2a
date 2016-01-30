class stuff(object):
    pass
pkg = stuff()
pkg.a = 'a'
pkg.b = 'b'
pkg.c = 'c'
pkg.d = 13
pkg.e = True

print vars(pkg)
print vars(pkg)['e']


for key,val in vars(pkg).items():
    exec(key + '=val')
print d
