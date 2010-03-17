
import sets
set = sets.Set

phi = set()

natural_numbers = [phi.copy()]
nats = natural_numbers

N = 20

for k in range(N):
    next = set()
    for j in range(len(nats)):
        next.add(nats[j].copy())
    nats.append(next)
    
for k in range(len(nats)):
    print len(nats[k])

for k in range(5):
    print len(nats[k]), "\t", nats[k]
    
