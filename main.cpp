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
			if(backbone[i][j])
				res.push_back({i, j});
	return res;
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

    printf("dimensions: %d x %d\n", n, m);
    int total_target = 0;
    for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			total_target += board[i][j] == '.';
    printf("total targets: %d\n", total_target);
    int area = (2 * radius + 1) * (2 * radius + 1);
    printf("routers to cover: %d\n", (total_target + area - 1) / area);

    printf("available routers: %d\n", budget / cost_router);

    printf("upper bound for score: %dM\n", ((long long)total_target * 1000LL + (long long)budget) / 1000000LL);

	for(auto i : make_backbone({{0,0}, {1,4}, {5,5}, {3,3}, {6,7}}))
		printf("%d %d\n", i.i, i.j);
	
    return 0;
}
