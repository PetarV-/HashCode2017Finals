#include <iostream>
#include <queue>
#include <vector>

using namespace std;

const int N = 1000;
int n, m, radius, cost_edge, cost_router, budget, start_i, start_j;
char board[N][N];

struct coord
{
    int i, j;
} ;

// ----------
// backbone

struct edge
{
	int from, to;
	int len;
} ;

bool operator<(edge a, edge b)
{
	// > for priority_queue as min-heap
	return a.len == b.len ? a.from == b.from ? a.to < b.to : a.from < b.from : a.len > b.len;
}

int dist(coord a, coord b)
{
	return max(abs(a.i - b.i), abs(a.j - b.j));
}

int sgn(int x) { return x ? (x > 0 ? 1 : -1) : 0; }

bool in_mst[N];
bool backbone[N][N];
vector<coord> make_backbone(vector<coord> routers)
{
    routers.push_back({start_i, start_j});
    
    priority_queue<edge> edges;
	for(int i = 0; i < routers.size(); i++)
		in_mst[i] = false;

	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			backbone[i][j] = false;
	
	in_mst[0] = true;
	for(int i = 1; i < routers.size(); i++)
		edges.push({0, i, dist(routers[0], routers[i])});

	for(int iter = 1; iter < routers.size(); iter++)
	{
		while(in_mst[edges.top().to])
			edges.pop();

		edge curr = edges.top();
		edges.pop();

		in_mst[curr.to] = true;
		for(int i = 0; i < routers.size(); i++)
			if(!in_mst[i])
				edges.push({curr.to, i, dist(routers[curr.to], routers[i])});

		// draw edge
		int ii = routers[curr.from].i, jj = routers[curr.from].j;
		backbone[ii][jj] = true;
		while(ii != routers[curr.to].i || jj != routers[curr.to].j)
		{
		    ii += sgn(routers[curr.to].i - ii);
		    jj += sgn(routers[curr.to].j - jj);
		    backbone[ii][jj] = true;
		}
	}
	
	vector<coord> res;
	for(int i = 0; i < n; i++)
	    for(int j = 0; j < n; j++)
		    if(backbone[i][j] && (i != start_i || j != start_j))
			    res.push_back({i, j});
	return res;
}

// ----------

int value(vector<coord> routers)
{
    bool mark[n][m];
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
        {
            mark[i][j] = false;
        }
    }

    int ret = 0;

    for (int i=0;i<routers.size();i++)
    {
        // check if outta bounds
        if (routers[i].i < 0 || routers[i].i > n - 1 || routers[i].j < 0 || routers[i].j > m - 1) return -2;
        // check if wallhacked
        if (board[routers[i].i][routers[i].j] == '#') return -1;
        // otherwise, expand up
        for (int dx=0;dx>=-radius;dx--)
        {
            int xt = routers[i].i + dx;
            if (xt < 0) break;
            
            // top-left
            for (int dy=0;dy>=-radius;dy--)
            {
                int yt = routers[i].j + dy;
                if (yt < 0) break;
                if (board[xt][yt] == '#') break; // no need to go on
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (board[xt][yt] == '.') ret++;
                }
            }

            // top-right
            for (int dy=0;dy<=radius;dy++)
            {
                int yt = routers[i].j + dy;
                if (yt > m - 1) break;
                if (board[xt][yt] == '#') break;
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (board[xt][yt] == '.') ret++;
                }
            }
        }

        // now expand down
        for (int dx=0;dx<=radius;dx++)
        {
            int xt = routers[i].i + dx;
            if (xt > n - 1) break;
            
            // bottom-left
            for (int dy=0;dy>=-radius;dy--)
            {
                int yt = routers[i].j + dy;
                if (yt < 0) break;
                if (board[xt][yt] == '#') break; // no need to go on
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (board[xt][yt] == '.') ret++;
                }
            }

            // top-right
            for (int dy=0;dy<=radius;dy++)
            {
                int yt = routers[i].j + dy;
                if (yt > m - 1) break;
                if (board[xt][yt] == '#') break;
                if (!mark[xt][yt])
                {
                    mark[xt][yt] = true;
                    if (board[xt][yt] == '.') ret++;
                }
            }
        }
    }
    return ret;
}

// -----------
// output

void write_output(string filename, vector<coord> routers, vector<coord> backbone)
{
    FILE *f = fopen(filename.c_str(), "w");
    if(!f)
    {
		fprintf(stderr, "Opening file %s for writing failed!\n", filename.c_str());
		return;
    }

    fprintf(f, "%d\n", backbone.size());
	for(coord i : backbone)
		fprintf(f, "%d %d\n", i.i, i.j);

	fprintf(f, "%d\n", routers.size());
	for(coord i : routers)
		fprintf(f, "%d %d\n", i.i, i.j);

	fclose(f);
}

// ----------

int main()
{
    scanf("%d %d %d", &n, &m, &radius);
    scanf("%d %d %d", &cost_edge, &cost_router, &budget);
    scanf("%d %d", &start_i, &start_j);

    for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			scanf(" %c", &board[i][j]);
	
    return 0;
}
