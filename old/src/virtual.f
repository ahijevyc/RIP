c                                                                     c
c*********************************************************************c
c                                                                     c
      function virtual(temp,ratmix)
c
c   This function returns virtual temperature in K, given temperature
c      in K and mixing ratio in kg/kg.
c
      include 'comconst'
c
      virtual=temp*(eps+ratmix)/(eps*(1.+ratmix))
      return
      end
