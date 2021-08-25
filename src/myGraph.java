import java.util.*;

public class myGraph
{
	protected String[] names;	// 1-d array to store the vertices
	protected boolean[][] Edges;	// 2-d array to store adjacencies between vertices,

	protected int numVertices;	
	protected int numEdges;

	// Default constructor. Sets aside enough capacity for one vertex
	public myGraph( )		
	{			
		this(1);
	}

	// Constructor that sets aside as much capacity as specified by the user
	public myGraph(int capacity)	
	{			
		names = new String[capacity];
		Edges = new boolean[capacity][];

		// We use only the portion of the matrix below the main diagonal to store the edges
		for(int i = 0; i < Edges.length; i++)
		{
			Edges[i] = new boolean[i+1];			
			for(int j = 0; j < i; j++)
				Edges[i][j] = false;
		}
	}
	
	public int numberOfVertices()
   	{
       	         return numVertices;
        }

        public int numberOfEdges()
        {
       	         return numEdges;
        }

	// Finds the location at which a vertex is stored in Vertices. 
	// Returns -1 if vertex not found
	public int getIndex(String vertex)
	{
		for(int i = 0; i < numVertices; i++)
			if(vertex.equals(names[i]))
				return i;

		return -1;
	}

	// Resizes the array of vertices. Can make it larger or smaller, depending
	// on what newSize is. 
	protected String[] resize(String[] array, int newSize)
	{
		String[] temp = new String[newSize];
		
		int smallerSize = newSize;
		if(array.length < smallerSize)
			smallerSize = array.length;

		for(int i = 0; i < smallerSize; i++)
			temp[i] = array[i];

		return temp;
	}
			


        // Resizes the array of adjacencies. Can make it larger or smaller, depending
        // on what newSize is.
        protected boolean[][] resize(boolean[][] array, int newSize)
        {
                boolean[][] temp = new boolean[newSize][];
        	int smallerSize = newSize;
        	if(array.length < smallerSize)
        		smallerSize = array.length;

		for(int i = 0; i < newSize; i++)
		{
			temp[i] = new boolean[i+1];
        		for(int j = 0; j < i+1; j++)
			{
            			if (i < smallerSize)
					temp[i][j] = array[i][j];
				else
					temp[i][j] = false;
			}
        	}
        	return temp;
    	}


	// Adds a new vertex
	public void addVertex(String newVertex)
	{
		if(getIndex(newVertex) != -1)
		{
			System.out.print("addVertex: ");
			System.out.print(newVertex);
			System.out.println(" failed -- vertex already exists.");
			return;
		}
		
		// if array of vertices is full, we need to expand it and 
		// also expand Edges
		if (names.length == numVertices)
		{
			names = resize(names, 2*numVertices+1);
			Edges = resize(Edges, 2*numVertices+1);
		}
			
		names[numVertices++] = newVertex;
	}


	// delete the given vertex
	public void deleteVertex(String vertex)
	{
		int i = getIndex(vertex);
		if(i == -1)
		{
			System.out.print("deleteVertex: ");
			System.out.print(vertex);
			System.out.println(" failed -- it does not exist.");
			return;
		}

		// Move the last vertex up to fill the hole made by vertex i
		names[i] = names[numVertices-1];
		names[numVertices-1] = null;

		for(int j = 0; j < i; j++)
		{
			// For every edge incident on vertex i, decrement numEdges
			if(Edges[i][j])
				numEdges--;
			// Move the first i elements of the last row in the adj. matrix into the ith row		
			Edges[i][j] = Edges[numVertices-1][j];
		}
		
		for(int j = i+1; j < numVertices-1; j++)
		{
			// For every edge incident on vertex i, decrement numEdges
			if(Edges[j][i])
				numEdges--;
			// Move the remaining elements of the last row in the adj. matrix into the ith collumn
			Edges[j][i] = Edges[numVertices-1][j];
		}
		
		// Delete the last row in the matrix
		for(int k = 0; k<numVertices; k++)
			Edges[numVertices-1][k] = false;

		numVertices--;

	}
	// Adds a new edge
	public void addEdge(String vertex1, String vertex2)
	{
		int i = getIndex(vertex1);
		if(i == -1)
		{
			System.out.print("addEdge failed: ");
			System.out.print(vertex1);
			System.out.println(" does not exist.");
            		return;
        	}

		int j = getIndex(vertex2);
		if(j == -1)
		{
			System.out.print("addEdge failed: ");
			System.out.print(vertex2);
			System.out.println(" does not exist.");
            		return;
        	}

        	if (i>=j) Edges[i][j] = true;
        	else Edges[j][i] = true;

		numEdges++;
	}


	// deletes a given edge
	public void deleteEdge(String vertex1, String vertex2)
	{
		int i = getIndex(vertex1);
        if(i == -1)
        {
            System.out.print("deleteEdge failed: ");
            System.out.print(vertex1);
            System.out.println(" does not exist.");
            return;
        }

        int j = getIndex(vertex2);
        if(j == -1)
        {
            System.out.print("deleteEdge failed: ");
            System.out.print(vertex2);
            System.out.println(" does not exist.");
            return;
        }


        if (i>=j) 
        {
        	if (Edges[i][j])
        	{
        		Edges[i][j] = false;
        		numEdges--;
        	}
        }
        else 
        {
        	if (Edges[j][i])
        	{
        		Edges[j][i] = false;
        		numEdges--;
        	}
        }
	}

	// returns the names of all the neighbors of a given vertex in an array of Strings
	public String[] getNeighbors(String vertex)
	{
		String[] neighbors = new String[numVertices];
		int numNeighbors = 0;
		
		int source = getIndex(vertex);
		if(source == -1)
		{
			System.out.print("getNeighbors failed: Vertex ");
            		System.out.print(vertex);
            		System.out.println(" does not exist.");
            		return new String[0];
        	}

		for(int j = 0; j < numVertices; j++)
		{
			boolean edge = false;
			if (j <= source) edge = Edges[source][j];
			else edge = Edges[j][source];
			
			if(edge)
				neighbors[numNeighbors++] = new String(names[j]);
		}
		
		neighbors = resize(neighbors, numNeighbors);
		return neighbors;
	}

	// returns the degree of a vertex with given name
    	public int degree(String vertex)
    	{
		// Get the index of the vertex
		int i = getIndex(vertex);
		if(i == -1)
		{
			System.out.print("In degree: ");
			System.out.print(vertex);
			System.out.println(" does not exist.");
			return -1;
		}

		// Call the other degree function that returns the degree
		// of a vertex, given its index
		return degree(i);
    	}

		
	// returns the degree of a vertex with given index
	public int degree(int index)
	{
        	int numNeighbors = 0;

		// Scan the row corresponding to vertex in the adjacency
		// matrix, counting the number of neighbors
        	for (int j = 0; j <= index; j++)
        		if(Edges[index][j])
            			numNeighbors++;
        
        	for (int j = index+1; j < numVertices; j++)
        		if(Edges[j][index])
        			numNeighbors++;

		return numNeighbors;	
	}

    	// returns the indices of all the neighbors of a given vertex with index
    	public int[] getNeighbors(int index)
    	{
        	int numNeighbors = degree(index);
        	int[] neighbors = new int[numNeighbors];

		// Scan the row corresponding to vertex in the adjacency matrix 
        	numNeighbors = 0;
        
        	for(int j = 0; j < numVertices; j++)
        	{
        		boolean edge = false;
        		if (j <= index) edge = Edges[index][j];
				else edge = Edges[j][index];

        		if(edge)
            			neighbors[numNeighbors++] = j;
        	}
        	return neighbors;
    	}

   	/* This is the depth first traversal function discussed in class.
	   It returns an integer array representing the depth first tree.
	   The function uses the Stack defined in the Java Collections.  
	   Data structures defined in Java Collections can be accessed by
	   importing java.util.* (as is done at the top of this file).
	   Data structures defined in Java Collections are generic in
	   the sense that one can store any Object into them. 
	   Notice the way a Stack is defined in the code below:
	 	   Stack<Integer> s = new Stack<Integer>();
	   Inside the angle brackets I am explicitly stating the kind of
	   objects I want stored in the stack. This is a new feature
	   of Java (introduced ini JDK 1.5) called generics. The advantage
	   of explicitly stating the type of objects that will go into the 
	   stack up-front is that when we do a pop or a peek on the stack
	   we get an Integer object back. Without the use of generics,
	   we would have gotten an object back and would have to typecast
	   it into an Integer (an action that could lead to run-time errors).
	   It is possible that at home you may be using an older compiler
	   that does not support generics. In this case this code will
	   not compile. To make your code compile, replace the above stack
	   definition by this: Stack s = new Stack(); Also, the line of code
	   that peeks into the stack needs to be changed to 
		int currentVertex = ((Integer)(s.peek())).intValue();
	   Notice the explicit typecasting of what is returned. The other
	   lines of code involving the Stack stay the same.		*/


        public int[] depthFirstTraversal(String source)
        {

                // Getting the index of the source vertex and
                // checking if the vertex really exists
                int sourceIndex = getIndex(source);
                if(sourceIndex == -1)
                {
                        System.out.print("In depthFirstTraversal: vertex ");
                        System.out.print(source);
                        System.out.println(" is missing.");
                        return null;
                }

                // Defining and initializing the visited array
                boolean[] visited = new boolean[numVertices];
                for(int j = 0; j < numVertices; j++)
                        visited[j] = false;

                // Defining and initializing the stack
                Stack<Integer> s = new Stack<Integer>();

                // Defining and initializing the depth first traversal tree
                int[] dFTTree = new int[numVertices];
                for(int j = 0; j < numVertices; j++)
                        dFTTree[j] = -1;

                boolean more;
                do
                {
			
			// Marking the source as visited and pushing the source
			// into the stack
                	visited[sourceIndex] = true;
                	s.push(sourceIndex);

                        // The traversal can go on while the stack
                        // contains a vertex to process
                        while(!s.empty())
                        {
                                // Peek at the current vertex
                                int currentVertex = (s.peek()).intValue();

                                System.out.println(names[currentVertex]);

                                // Get the indices of the neighbors of the current vertex
                                int[] neighbors = getNeighbors(currentVertex);

                                // Scan the neighbors of the current vertex, looking
                                // for an unvisited neighbor
                                int j = 0;
                                boolean found = false;
                                while (j < neighbors.length && !found)
                                {
                                        // If an unvisited neighbor has been found,
                                        // then push it into the stack, mark it visited,
                                        // make the current vertex its parent,
                                        // and get out of the while-loop by setting found
                                        // to true
                                        if(!visited[neighbors[j]])
                                        {
                                                s.push(neighbors[j]);
                                                visited[neighbors[j]] = true;
                                                found = true;
                                                dFTTree[neighbors[j]] = currentVertex;
                                        }

                                        j++; // scan the next vertex
                                }

                                // If no unvisited vertices have been found, it is time
                                // to backtrack
                                if (!found)
                                        s.pop();
                        } // end of while-stack-not-empty loop

                        // Determine if there are more unvisited vertices
                        // by scanning the visited array and looking for the
                        // first unvisited vertex
                        more = false;
                        int j = 0;
                        while(j < numVertices && !more)
                        {
                                if(!visited[j])
				{
                                        more = true;
					sourceIndex = j;
				}
                                j++;
                        }

                }
                while(more);

                return dFTTree;

        } // end of function

	public int[] recursiveDepthFirstTraversal(int currentVertex, boolean [] visited, int [] dFTTree)
	{
		System.out.println(names[currentVertex]);

		// Get the indices of the neighbors of the current vertex
		int[] neighbors = getNeighbors(currentVertex);

		// Scan the neighbors of the current vertex, looking
		// for an unvisited neighbor
		int j = 0;
		while (j < neighbors.length)
		{
			// If an unvisited neighbor has been found,
			// then push it into the stack, mark it visited,
			// make the current vertex its parent,
			// and get out of the while-loop by setting found
			// to true
			if(!visited[neighbors[j]])
			{
				visited[neighbors[j]] = true;
				dFTTree[neighbors[j]] = currentVertex;
				recursiveDepthFirstTraversal(neighbors[j], visited, dFTTree);

				// Re-output the current vertex for tracing purposes
				System.out.println(names[currentVertex]);
			}

			j++; // scan the next vertex
		}

		// We've visited all our children, so return the tree.
		return dFTTree;
	}

        public int[] rDepthFirstTraversal(String source)
        {
                // Getting the index of the source vertex and
                // checking if the vertex really exists
                int sourceIndex = getIndex(source);
                if(sourceIndex == -1)
                {
                        System.out.print("In rDepthFirstTraversal: vertex ");
                        System.out.print(source);
                        System.out.println(" is missing.");
                        return null;
                }

                // Defining and initializing the visited array
                boolean[] visited = new boolean[numVertices];
                for(int j = 0; j < numVertices; j++)
                        visited[j] = false;

                // Defining and initializing the depth first traversal tree
                int[] dFTTree = new int[numVertices];
                for(int j = 0; j < numVertices; j++)
                        dFTTree[j] = -1;

                boolean more;
                do
                {
	                visited[sourceIndex] = true;
			recursiveDepthFirstTraversal(sourceIndex, visited, dFTTree);

                        // Determine if there are more unvisited vertices
                        // by scanning the visited array and looking for the
                        // first unvisited vertex
                        more = false;
                        int j = 0;
                        while(j < numVertices && !more)
                        {
                                if(!visited[j])
				{
                                        more = true;
					sourceIndex = j;
				}
                                j++;
                        }

                }
                while(more);

                return dFTTree;

        } // end of function

} // end of class

