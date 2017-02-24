def cs_prop(nodes, elem):
    # Calculate cross sectional properties
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    nele = len(elem[0])
    node = elem[0]+elem[1]
    nnode = 0
    j = 0
    
    while node:
        i = [ii for ii,x in enumerate(node) if x==node[0]]
        for ii in sorted(i, reverse=True):
            del node[ii]
        if len(i)==2:
            j += 1
        nnode += 1
    
    # classify the section type
    if j == nele:
        section = 'close' #single cell
    elif j == nele-1:
        section = 'open' #singly-branched
    else:
        #connected
        section = 'open' #multi-branched
        #disconnected
        #in the future it would be good to handle multiple cross-sections in
        #one model, etc., for now the code will bomb if more than a single
        #section is used - due to inability to calculate section properties.
        #section = 'arbitrary'; #arbitrary section unidentified
    
    # if the section is close re-order the element 
    # The code is not translated to python, to be completed
    # The code works if the elements are given in a sequential order
    # In order to include randomly given elements, the sorting script needs to be writen
    #if section=='close':
    #    xnele = (nele-1)
    #    for i in range(xnele)
    #        en = elem
    #        en[i][1] = 0
    #        [m,n] = find(elem(i,2)==en(:,1:2))
    #        if n==1
    #            elem(i+1,:) = en(m,:)
    #            elem(m,:) = en(i+1,:)
    #        elseif n == 2
    #            elem(i+1,:) = en(m,[2 1 3])
    #            elem(m,:) = en(i+1,[2 1 3])
    
    # find the element properties
    tt = []
    xm = []
    ym = []
    xd = []
    yd = []
    L = []
    for i in range(nele):
        sn = elem[0][i]
        fn = elem[1][i]
        # thickness of the element
        tt = tt+[elem[2][i]]
        # compute the coordinate of the mid point of the element
        xm = xm + [mean([nodes[0][sn], nodes[0][fn]])]
        ym = ym + [mean([nodes[1][sn], nodes[1][fn]])]
        # compute the dimension of the element
        xd = xd + [(nodes[0][fn]-nodes[0][sn])]
        yd = yd + [(nodes[1][fn]-nodes[1][sn])]
        # compute the length of the element
        L = L + [sqrt(xd[i]**2+yd[i]**2)]
    
    # compute the cross section area
    A = sum([x*t for x in L])
    # compute the centroid 
    xc = sum([a*b*c for a,b,c in zip(L,tt,xm)])/A
    yc = sum([a*b*c for a,b,c in zip(L,tt,ym)])/A
    
    if abs(xc/sqrt(A)) < 1e-12:
        xc = 0
    
    if abs(yc/sqrt(A)) < 1e-12:
        yc = 0
    
    # compute the moment of inertia
    Ix = sum([sum(a) for a in zip([a**2*b*c/12 for a,b,c in zip(yd,L,tt)], [(a-yc)**2*b*c for a,b,c in zip(ym,L,tt)])])
    Iy = sum([sum(a) for a in zip([a**2*b*c/12 for a,b,c in zip(xd,L,tt)], [(a-xc)**2*b*c for a,b,c in zip(xm,L,tt)])])
    Ixy = sum([sum(a) for a in zip([a*b*c*d/12 for a,b,c,d in zip(xd,yd,L,tt)], [(a-xc)*(b-yc)*c*d for a,b,c,d in zip(xm,ym,L,tt)])])
    
    if abs(Ixy/A**2) < 1e-12:
        Ixy = 0
    
    # compute the rotation angle for the principal axes
    theta_principal = atan((-2*Ixy)/(Ix-Iy))/2
    
    # transfer the section coordinates to the centroid principal coordinates
    coord12 = [[a-xc for a in nodes[0]],[a-yc for a in nodes[1]]]
    coord12 = np.array([[cos(theta_principal), sin(theta_principal)],[-sin(theta_principal), cos(theta_principal)]]).dot(nodes)
    
    # find the element properties
    for i in range(nele):
        sn = elem[0][i]
        fn = elem[1][i]
        # compute the coordinate of the mid point of the element
        xm = xm + [mean([coord12[0][sn], coord12[0][fn]])]
        ym = ym + [mean([coord12[1][sn], coord12[1][fn]])]
        # compute the dimension of the element
        xd = xd + [(coord12[0][fn]-coord12[0][sn])]
        yd = yd + [(coord12[1][fn]-coord12[1][sn])]
    
    # compute the principal moment of inertia
    I1 = sum([sum(a) for a in zip([a**2*b*c/12 for a,b,c in zip(yd,L,tt)], [(a-yc)**2*b*c for a,b,c in zip(ym,L,tt)])])
    I2 = sum([sum(a) for a in zip([a**2*b*c/12 for a,b,c in zip(xd,L,tt)], [(a-xc)**2*b*c for a,b,c in zip(xm,L,tt)])])
    return I1