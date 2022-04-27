#pascal triangle
def sortNumbers (Finaline):
    Valeur = Finaline
    if Valeur [0] > Valeur [1]:
        Valeur [0], Valeur [1] = Valeur [1], Valeur [0]
    if Valeur [1] > Valeur [2]:
        Valeur [1], Valeur [2] = Valeur [2], Valeur [1]
    return Valeur

def find_three_largest (Finaline):
    if len(sortNumbers) < 3:
        return None
    largest = sortNumbers(Finaline[:3])
    if len(sortNumbers) > 3 :
        return largest
    
    for n in Finaline [3:] :
        if largest[0] < n :
            largest[0] = n
        largest = sortNumbers(Finaline) 
    return largest