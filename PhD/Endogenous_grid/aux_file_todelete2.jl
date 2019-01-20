function rouwenhorst(;
                    rho = 0.9,
                    p = (1.0+0.9)/2, #p = (1.0+rho)/2,
                    q = p,
                    sigma_e = 0.1, # Dont know
                    sigma_y = sqrt((sigma_e^2)/(1-rho^2)),
                    ny = 4, # NUmber of states
                    my = 0) # Who knows) # Number of states

        phy = sigma_y*sqrt(ny-1)
        ymax = phy
        ymin = -phy
        ygrid = range(ymin, ymax, legnth=ny)

        Piy2 = [[p 1.0-p];[1.0-q q]]
        Piyn1 = copy(Piy2)

        for jj = 1:(ny-2)
            num_rows = size(Piyn1,1)
            mat1     = zeros(num_rows+1, num_rows+1)
            mat2, mat3, mat4 = copy(mat1), copy(mat1), copy(mat1)

            mat1[1:end-1, 1:end-1]  = Piyn1
            mat2[1:end-1, 2:end]    = Piyn1
            mat3[2:end, 1:end-1]    = Piyn1
            mat4[2:end, 2:end]      = Piyn1

            Piyn1 = p*mat1 + (1-p)*mat2 + (1-q)*mat3 + q*mat4
            Pyn1[2:end-1, :] = Piyn1[2:end-1, :] / 2
        end
        Piy     = copy(Piyn1)
        Piy_aux = copy(Piy')
        vals = eigvals(Piy_aux)
        vecs = eigvecs(Piy_aux)
        todelete, ind_val  = findmin(abs.(vals.-1.0))
        Pyinv       = vecs[:, ind_val]/sum(vecs[:, ind_val])

        sum(Pyinv.>=0.0) == ny || throw(error("Negative elements in invariant distribution"))
        egrid   = exp.(ygrid + my*ones(ny))
        #return (Pyinv=Pyinv, Piy=Piy ,egrid=egrid ,ny=ny)
        return Pyinv
           # Persistence of the income process
end
