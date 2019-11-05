
#include "bounding_sphere.hpp"
#include <stack>

template<typename HullType> void
KDTree<HullType>::construct()
{
    m_nodes.clear();
    m_hulls.clear();
    if (m_lst.empty()) return;

    std::iota(m_lst.begin(), m_lst.end(), 0);

    // Determine bounding box of considered domain.
    auto box = Eigen::AlignedBox3d{};
    for (auto i = 0u; i < m_lst.size(); ++i)
        box.extend(entityPosition(i));

    auto ni = addNode(0, static_cast<unsigned int>(m_lst.size()));
    construct(ni, box, 0, static_cast<unsigned int>(m_lst.size()));
}

template<typename HullType> void
KDTree<HullType>::construct(unsigned int node, Eigen::AlignedBox3d const& box, unsigned int b,
    unsigned int n)
{
    // If only one element is left end recursion.
    //if (n == 1) return;
    if (n < 10) return;

    // Determine longest side of bounding box.
    auto max_dir = 0;
    auto d = box.diagonal().eval();
    if (d(1) >= d(0) && d(1) >= d(2))
        max_dir = 1;
    else if (d(2) >= d(0) && d(2) >= d(1))
        max_dir = 2;

#ifdef _DEBUG
    for (auto i = 0u; i < n; ++i)
    {
        if (!box.contains(entityPosition(m_lst[b + i])))
            std::cerr << "ERROR: Bounding box wrong!" << std::endl;
    }
#endif

    // Sort range according to center of the longest side.
    std::sort(m_lst.begin() + b, m_lst.begin() + b + n, 
        [&](unsigned int a, unsigned int b)
        {
            return entityPosition(a)(max_dir) < entityPosition(b)(max_dir);
        }
    );

    auto hal = n / 2;
    auto n0 = addNode(b      , hal    );
    auto n1 = addNode(b + hal, n - hal);
    m_nodes[node].children[0] = n0;
    m_nodes[node].children[1] = n1;

    auto c = 0.5 * (
        entityPosition(m_lst[b + hal -1])(max_dir) + 
        entityPosition(m_lst[b + hal   ])(max_dir));
    auto l_box = box; l_box.max()(max_dir) = c;
    auto r_box = box; r_box.min()(max_dir) = c;

    construct(m_nodes[node].children[0], l_box, b, hal);
    construct(m_nodes[node].children[1], r_box, b + hal, n - hal);
}

template<typename HullType> void
KDTree<HullType>::traverseDepthFirst(TraversalPredicate pred, TraversalCallback cb,
    TraversalPriorityLess const& pless) const
{
    if (m_nodes.empty())
        return;

    if (pred(0, 0)) 
        traverseDepthFirst(0, 0, pred, cb, pless);
}

template<typename HullType> void
KDTree<HullType>::traverseDepthFirst(unsigned int node_index, 
    unsigned int depth, TraversalPredicate pred, TraversalCallback cb,
    TraversalPriorityLess const& pless) const
{

    //auto pending = std::stack<QueueItem>{};
    //pending.push({node_index, depth});
    //while (!pending.empty())
    //{
    //    auto n = pending.top().n;
    //    auto d = pending.top().d;
    //    auto const& node = m_nodes[n];
    //    pending.pop();

    //    cb(n, d);
    //    auto is_pred = pred(n, d);
    //    if (!node.is_leaf() && is_pred)
    //    {
    //        if (pless && !pless(node.children))
    //        {
    //            pending.push({ static_cast<unsigned int>(node.children[1]), d + 1 });
    //            pending.push({ static_cast<unsigned int>(node.children[0]), d + 1 });
    //        }
    //        else
    //        {
    //            pending.push({ static_cast<unsigned int>(node.children[0]), d + 1 });
    //            pending.push({ static_cast<unsigned int>(node.children[1]), d + 1 });
    //        }
    //    }
    //}


    Node const& node = m_nodes[node_index];

    cb(node_index, depth);
    auto is_pred = pred(node_index, depth);
    if (!node.isLeaf() && is_pred)
    {
        if (pless && !pless(node.children))
        {
            traverseDepthFirst(m_nodes[node_index].children[1], depth + 1, pred, cb, pless);
            traverseDepthFirst(m_nodes[node_index].children[0], depth + 1, pred, cb, pless);
        }
        else
        {
            traverseDepthFirst(m_nodes[node_index].children[0], depth + 1, pred, cb, pless);
            traverseDepthFirst(m_nodes[node_index].children[1], depth + 1, pred, cb, pless);
        }
    }




        //    auto n = pending.front().n;
        //auto d = pending.front().d;
        //auto const& node = m_nodes[n];
        //pending.pop();

        //cb(n, d);
        //auto is_pred = pred(n, d);
        //if (!node.is_leaf() && is_pred)
        //{
        //    if (pless && !pless(node.children))
        //    {
        //        pending.push({ static_cast<unsigned int>(node.children[1]), d + 1 });
        //        pending.push({ static_cast<unsigned int>(node.children[0]), d + 1 });
        //    }
        //    else
        //    {
        //        pending.push({ static_cast<unsigned int>(node.children[0]), d + 1 });
        //        pending.push({ static_cast<unsigned int>(node.children[1]), d + 1 });
        //    }
        //}
}


template <typename HullType> void
KDTree<HullType>::traverseBreadthFirst(TraversalPredicate const& pred, 
    TraversalCallback const& cb, unsigned int start_node, TraversalPriorityLess const& pless,
    TraversalQueue& pending) const
{
    //auto pending = TraversalQueue{};

    cb(start_node, 0);
    if (pred(start_node, 0)) pending.push({ start_node, 0 });
    traverseBreadthFirst(pending, pred, cb, pless);
}

template <typename HullType> unsigned int
KDTree<HullType>::addNode(unsigned int b, unsigned int n)
{
    HullType hull;
    computeHull(b, n, hull);
    m_hulls.push_back(hull);
    m_nodes.push_back({ b, n });
    return static_cast<unsigned int>(m_nodes.size() - 1);
}

template <typename HullType> void
KDTree<HullType>::update()
{
    traverseDepthFirst(
        [&](unsigned int, unsigned int) { return true; },
        [&](unsigned int node_index, unsigned int)
        {
            auto const& nd = node(node_index);
            computeHull(nd.begin, nd.n, hull(node_index));
        }
    );
}

template <typename HullType> void 
KDTree<HullType>::traverseBreadthFirst(TraversalQueue& pending, 
    TraversalPredicate const& pred, TraversalCallback const& cb, TraversalPriorityLess const& pless) const
{
    while (!pending.empty())
    {
        auto n = pending.front().n;
        auto d = pending.front().d;
        auto const& node = m_nodes[n];
        pending.pop();

        cb(n, d);
        auto is_pred = pred(n, d);
        if (!node.is_leaf() && is_pred)
        {
            if (pless && !pless(node.children))
            {
                pending.push({ static_cast<unsigned int>(node.children[1]), d + 1 });
                pending.push({ static_cast<unsigned int>(node.children[0]), d + 1 });
            }
            else
            {
                pending.push({ static_cast<unsigned int>(node.children[0]), d + 1 });
                pending.push({ static_cast<unsigned int>(node.children[1]), d + 1 });
            }
        }
    }
}
