  {
      int grpctr, subctr;
      
      printf ( " number of groups: %d\n",top.atoms.ngrpname);
      
      for (grpctr=0; grpctr<top.atoms.ngrpname ; grpctr++) {
	  
	  printf ( "\t %3d subgroups: %2d \n", grpctr, top.atoms.grps[grpctr].nr);
	  for ( subctr=0; subctr<top.atoms.grps[grpctr].nr && subctr<=10; subctr++) {
	      if ( ! top.atoms.grpname[grpctr][subctr]) break;
	      printf ( "\t\t %3d %2d  %-s \n", grpctr, subctr, top.atoms.grpname[grpctr][subctr]);
	  }
      }
  }
  exit (0);
