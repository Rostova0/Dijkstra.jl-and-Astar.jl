using LinearAlgebra
using DataStructures #to use PriorityQueue
using ImageView #to use imshow
using Colors
using Gtk

function toMatrix(namefile) #converts the data of a map into a matrix of char
    open("./"*namefile,"r") do file

        readline(file)
        
        hg::String= split(readline(file)," ")[2]
        wd::String= split(readline(file)," ")[2]
        height=parse(Int, hg)
        width=parse(Int, wd)
        readline(file)

        ma_matrix::Matrix{Char}=Matrix{Char}(undef,height,width)

        for i in 1:height
            lines=readline(file)
            for j in 1:width
                ma_matrix[i,j]=lines[j]
            end
        end   
        return ma_matrix
    end    
end

cost::Dict=Dict{Char,Int64}( #assigns a cost to each area of the map
            'S'=>5,
            'W'=>8,
            'G'=>1,
            '.'=>1,
            'T'=>99999,
            '@'=>99999,
            'O'=>99999)

function inMatrix(point::Tuple{Int64,Int64},h::Int64,w::Int64) #test if a point is in the matrix
    if point[1]>0 && point[2]>0
        if point[1]<=h && point[2]<=w
            return true
        end
    end
    return false
end

function heuristic(xd::Int64,yd::Int64,xa::Int64,ya::Int64) #heuristic Manhattan distance
    return abs(xd-xa)+abs(yd-ya)
end

function Astar(map,height::Int64,width::Int64,visited::Matrix{Bool},dist::Matrix{Int64},prio::PriorityQueue{Tuple{Int64,Int64},Int64},pred::Matrix{Tuple{Int64,Int64}},A::Tuple{Int64,Int64}) # Function that implements Astar
    
    nb_visit::Int64=1 #To get the number of states evaluated
    matrx=toMatrix(map)
    while visited[A[1],A[2]]!=true && length(prio)!=0
        p=dequeue!(prio) #picks the one with the lowest priority
        x=p[1]
        y=p[2]
        visited[x,y]=true
        nb_visit+=1
        # Update dist,pred and priority value of the adjacent state
        # of the picked state only if the current
        # distance is greater than new distance
        if inMatrix((x,y+1),height,width) && visited[x,y+1]==false && cost[matrx[x,y+1]]!=99999
            dn=dist[x,y]+cost[matrx[x,y+1]]
            if dn<dist[x,y+1]
                prio[(x,y+1)]=dn+heuristic(x,(y+1),A[1],A[2]) #Change an existing item’s priority to a lower number adding the heuristic
                dist[x,y+1]=dn
                pred[x,y+1]=p
            end
        end

        if inMatrix((x-1,y),height,width) && visited[x-1,y]==false && cost[matrx[x-1,y]]!=99999
            dw=dist[x,y]+cost[matrx[x-1,y]]
            if dw<dist[x-1,y]
                prio[(x-1,y)]=dw+heuristic((x-1),y,A[1],A[2])
                dist[x-1,y]=dw
                pred[x-1,y]=p
            end
        end

        if inMatrix((x+1,y),height,width) && visited[x+1,y]==false && cost[matrx[x+1,y]]!=99999
            de=dist[x,y]+cost[matrx[x+1,y]]
            if de<dist[x+1,y]
                prio[(x+1,y)]=de+heuristic((x+1),y,A[1],A[2])
                dist[x+1,y]=de
                pred[x+1,y]=p

            end
        end

        if inMatrix((x,y-1),height,width) && visited[x,y-1]==false && cost[matrx[x,y-1]]!=99999
            ds=dist[x,y]+cost[matrx[x,y-1]]
            if ds<dist[x,y-1]
                prio[(x,y-1)]=ds+heuristic(x,(y-1),A[1],A[2])
                dist[x,y-1]=ds
                pred[x,y-1]=p
            end
        end
    end
    return nb_visit
end

function reconstruct_path(prev::Matrix{Tuple{Int64,Int64}},D::Tuple{Int64,Int64},A::Tuple{Int64,Int64})
    path=Vector{Tuple{Int64,Int64}}()
    current::Tuple{Int64,Int64}=A

    if A[1]>0 && A[2]>0 && A[1]<=size(prev,1) && A[2]<=size(prev,2)

        while current!=D #stop when we reach the start point
            push!(path,current)
            current=prev[current[1],current[2]]
        end
    end
    return path #return the number of states evaluated
end

#assigns a color to a character
rgb::Dict{Char,RGB{Float64}} = Dict('S'=>RGB(240/255,230/255,140/255),
                                        'W'=>RGB(51/255,51/255,255/255),
                                        'G'=>RGB(158/255,158/255,164/255),
                                        '.'=>RGB(158/255,158/255,164/255),
                                        'T'=>RGB(11/255,135/255,41/255),
                                        '@'=>RGB(.0/255,3/255,1/255),
                                        'O'=>RGB(.0/255,3/255,1/255))

#converts a character matrix into an RGB matrix of the same size                                    
function toRGB(matrx::Matrix{Char},matrpd::Matrix{Bool},pred::Vector{Tuple{Int64,Int64}}) #build the path
    matrxrgb::Matrix{RGB{Float64}}=Matrix{RGB{Float64}}(undef,(size(matrx,1)),(size(matrx,2)))
    for i in 1:size(matrx,1)
        for j in 1:size(matrx,2)
            #all the evaluated states take the red color
            if matrpd[i,j]==true 
                matrxrgb[i,j]=RGB(244/255,20/255,13/255)
            else
                matrxrgb[i,j]=rgb[matrx[i,j]]
            end
       end  
    end 
    #draw the shortest path in yellow
    for j in 1:size(pred,1)
        matrxrgb[pred[j][1],pred[j][2]]=RGB(223/255,244/255,13/255)
    end
    return matrxrgb
end

function algoAstar(map,D::Tuple{Int64,Int64},A::Tuple{Int64,Int64}) #main function

    #initialization of variables
    matrx::Matrix{Char}=toMatrix(map)
    height=size(matrx,1)
    width=size(matrx,2)
    visited::Matrix{Bool}=fill(false,height,width) #sets all elements of the matrix to false
    dist::Matrix{Int64}=fill(99999,height,width)
    dist[D[1],D[2]]=0
    visited[D[1],D[2]]=true
    prio=PriorityQueue{Tuple{Int64,Int64},Int64}() #create an empty priotyqueue
    prio[(D[1],D[2])]=0 #adds the start state to the queue with priority 0
    pred::Matrix{Tuple{Int64,Int64}}=fill((1,1),height,width)

    println("Distance D → A  : ",Astar(map,height,width,visited,dist,prio,pred,A))
    println("Number of states evaluated: ",dist[A[1],A[2]])
    prec=reconstruct_path(pred,D,A)
    img=toRGB(matrx,visited,prec)
    win=imshow(img); #to view the created image
    resize!(win["gui"]["window"], 800, 800)
end
