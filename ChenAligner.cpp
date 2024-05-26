#include "ChenAligner.h"

void SegmentCouple::insert_noOverlap(unsigned int segl, unsigned int segr)
{
	if (m_begin == nullptr)
	{
		m_begin = new SCNode(segl, segr);
		m_coverlength = segr - segl + 1;
		minl = segl;
		maxr = segr;
	}
	else
	{
		minl = std::min(minl, segl);
		maxr = std::max(maxr, segr);
		SCNode* ptr = m_begin;
		SCNode* lastp = nullptr;

		while (true)
		{
			if (ptr == nullptr)
			{
				if (lastp == nullptr)
				{
					m_begin = new SCNode(segl, segr);
					m_coverlength += segr - segl + 1;
				}
				else
				{
					lastp->next = new SCNode(segl, segr);
					m_coverlength += segr - std::max(segl, lastp->l) + 1;
				}
				return;
			}
			else if (segl > ptr->l)
			{
				if (segr > ptr->r)
				{
					lastp = ptr;
					ptr = ptr->next;
				}
				else
				{
					return;
				}
			}
			else if (segl == ptr->l)
			{
				if (segr > ptr->r)
				{
					if (lastp == nullptr)
						m_begin = ptr->next;
					else
						lastp->next = ptr->next;
					m_coverlength -= (ptr->next == nullptr ? ptr->r : std::min(ptr->r, ptr->next->l)) - (lastp == nullptr ? ptr->l : std::max(ptr->l, lastp->r)) + 1;
					SCNode* delptr = ptr;
					ptr = ptr->next;
					delete delptr;
				}
				else
				{
					return;
				}
			}
			else
			{
				if (segr >= ptr->r)
				{
					if (lastp == nullptr)
						m_begin = ptr->next;
					else
						lastp->next = ptr->next;
					m_coverlength -= (ptr->next == nullptr ? ptr->r : std::min(ptr->r, ptr->next->l)) - (lastp == nullptr ? ptr->l : std::max(ptr->l, lastp->r)) + 1;
					SCNode* delptr = ptr;
					ptr = ptr->next;
					delete delptr;
				}
				else
				{
					if (lastp == nullptr)
					{
						m_begin = new SCNode(segl, segr);
						m_begin->next = ptr;
						m_coverlength += std::min(ptr->l, segr) - segl + 1;
					}
					else
					{
						lastp->next = new SCNode(segl, segr);
						lastp->next->next = ptr;
						m_coverlength += std::min(segr, ptr->l) - std::max(segl, ptr->r) + 1;
					}
					return;
				}
}
		}
	}
}

void SegmentCouple::insert(unsigned int segl, unsigned int segr)
{
	if (m_begin == nullptr)
	{
		m_begin = new SCNode(segl, segr);
		m_coverlength = segr - segl + 1;
		minl = segl;
		maxr = segr;
	}
	else
	{
		minl = std::min(minl, segl);
		maxr = std::max(maxr, segr);
		SCNode* ptr = m_begin;
		SCNode* lastp = nullptr;

		while (ptr != nullptr)
		{
			if (segl > ptr->r + 1)	// l on the right
			{
				lastp = ptr;
				ptr = ptr->next;
			}
			else if (segl < ptr->l)	// l on the left
			{
				if (segr < ptr->l - 1)
				{
					m_coverlength += segr - segl + 1;
					if (lastp != nullptr)
					{
						lastp->next = new SCNode(segl, segr);
						lastp->next->next = ptr;
					}
					else
					{
						m_begin = new SCNode(segl, segr);
						m_begin->next = ptr;
					}
				}
				else
				{
					m_coverlength += ptr->l - segl;
					ptr->l = segl;
					if (segr > ptr->r)
					{
						SCNode* nextp = ptr->next;
						while (nextp != nullptr && segr > nextp->r)
						{
							m_coverlength -= nextp->r - nextp->l + 1;
							nextp = nextp->next;
							delete ptr->next;
							ptr->next = nextp;
						}
						if (nextp == nullptr || segr < nextp->l - 1)
						{
							m_coverlength += segr - ptr->r;
							ptr->r = segr;
						}
						else
						{
							m_coverlength += nextp->l - ptr->r - 1;
							ptr->r = nextp->r;
							ptr->next = nextp->next;
							delete nextp;
						}
					}
				}
				break;
			}
			else	// l in the middle
			{
				if (segr > ptr->r)
				{
					SCNode* nextp = ptr->next;
					while (nextp != nullptr && segr > nextp->r)
					{
						m_coverlength -= nextp->r - nextp->l + 1;
						nextp = nextp->next;
						delete ptr->next;
						ptr->next = nextp;
					}
					if (nextp == nullptr || segr < nextp->l - 1)
					{
						m_coverlength += segr - ptr->r;
						ptr->r = segr;
					}
					else
					{
						m_coverlength += nextp->l - ptr->r - 1;
						ptr->r = nextp->r;
						ptr->next = nextp->next;
						delete nextp;
					}
				}
				break;
			}
		}

		if (ptr == nullptr)
		{
			m_coverlength += segr - segl + 1;
			lastp->next = new SCNode(segl, segr);

		}
	}

}

void SegmentCouple::print_noOverlap()
{
	std::cout << "Cover:" << m_coverlength << " Begin:" << minl << std::endl;

	SCNode* ptr = m_begin;
	while (ptr != nullptr)
	{
		for (int i = minl; i < ptr->l; i++)
		{
			std::cout << "_";
		}
		for (int i = ptr->l; i <= ptr->r; i++)
		{
			std::cout << "*";
		}
		for (int i = ptr->r + 1; i <= maxr; i++)
		{
			std::cout << "_";
		}
		std::cout << std::endl;
		ptr = ptr->next;
	}
}

void SegmentCouple::print()
{
	std::cout << "Cover:" << m_coverlength << " [" << minl << ", " << maxr << "]" << std::endl;

	bool do_print = true;
	SCNode* ptr = m_begin;
	for (unsigned int i = minl; i <= maxr; i++)
	{
		if (do_print)
		{
			std::cout << "*";
			if (i == ptr->r)
			{
				do_print = false;
				ptr = ptr->next;
			}
		}
		else
		{
			std::cout << "_";
			if (i == ptr->l - 1)
			{
				do_print = true;
			}
		}
	}
}

std::vector<unsigned int> SegmentCouple::best_joined()
{
	std::vector<unsigned int> dp;
	std::vector<SCNode*> hist;
	dp.resize(maxr - minl + 2, 0);
	hist.resize(maxr - minl + 2, nullptr);
	SCNode* ptr = m_begin;
	while (ptr != nullptr)
	{
		unsigned int maxlenwith = dp[ptr->l - minl] + ptr->r - ptr->l + 1;
		for (unsigned int i = ptr->r - minl + 1; i <= maxr - minl + 1; i++)
		{
			if (maxlenwith > dp[i])
			{
				dp[i] = maxlenwith;
				hist[i] = ptr;
			}
		}
		ptr = ptr->next;
	}

	std::vector<unsigned int> result;
	ptr = hist[maxr - minl + 1];
	while (ptr != nullptr)
	{
		result.push_back(ptr->r);
		result.push_back(ptr->l);
		ptr = hist[ptr->l - minl];
	}
	std::reverse(result.begin(), result.end());
	return result;
}

void SegmentCouple::best_joined(std::ostream& os)
{
	std::vector<unsigned int> res = best_joined();
	bool isBetween = true;
	for (auto it = res.begin(); it != res.end(); it++)
	{
		os << *it;
		if (isBetween)
			os << ":";
		else
			os << " ";
		isBetween = !isBetween;
	}
}

std::vector<std::pair<unsigned int, unsigned int>> SegmentCouple::get_all_segments() const
{
	SCNode* ptr = m_begin;
	std::vector<std::pair<unsigned int, unsigned int>> ans;
	while (ptr != nullptr)
	{
		ans.push_back(std::make_pair(ptr->l, ptr->r));
		ptr = ptr->next;
	}
	return ans;
}

std::vector<bool> SegmentCouple::get_individual_cover(unsigned int l, unsigned int r) const
{
	std::vector<bool> ans;
	ans.resize(r - l + 1, false);

	for (auto& seg : get_all_segments())
	{
		for (unsigned int i = seg.first; i <= seg.second; i++)
		{
			ans[i - l] = true;
		}
	}

	return ans;
}

void SegmentCouples::insert(unsigned int beginid, unsigned int endid, unsigned int segl, unsigned int segr)
{
	for (unsigned int i = beginid; i <= endid; i++)
	{
		auto it = storage.find(i);
		if (it == storage.end())
		{
			SegmentCouple sc;
			sc.insert(segl, segr);
			storage[i] = sc;
			max_coverlength = std::max(max_coverlength, storage[i].coverlength());
		}
		else
		{
			it->second.insert(segl, segr);
			max_coverlength = std::max(max_coverlength, it->second.coverlength());
		}
	}
}

std::vector<std::pair<unsigned int, SegmentCouple>> SegmentCouples::show_matches(unsigned int density_diff, std::vector<unsigned int> id_convertor)
{
	good_matches.clear();
	if (id_convertor.empty())
	{
		for (auto& sc : storage)
		{
			if (sc.second.coverlength() + density_diff >= max_coverlength)
			{
				good_matches.push_back(sc);

				if constexpr (confidence_requirement > 0)
				{
					if (good_matches.size() > confidence_requirement)
					{
						good_matches.clear();
						return good_matches;
					}
				}
			}
		}
	}
	else
	{
		for (auto& sc : storage)
		{
			if (sc.second.coverlength() + density_diff >= max_coverlength)
			{
				good_matches.push_back(std::make_pair(id_convertor[sc.first], sc.second));
			}
		}
	}
	std::sort(good_matches.begin(), good_matches.end(),
		[](const std::pair<unsigned int, SegmentCouple>& a, const std::pair<unsigned int, SegmentCouple>& b) {
			return a.second.coverlength() > b.second.coverlength();
		});

	return good_matches;
}

void SegmentCouples::insert(std::vector<unsigned int> id, unsigned int segl, unsigned int segr)
{
	for (unsigned int& i : id)
	{
		auto it = storage.find(i);
		if (it == storage.end())
		{
			SegmentCouple sc;
			sc.insert(segl, segr);
			storage[i] = sc;
			max_coverlength = std::max(max_coverlength, storage[i].coverlength());
		}
		else
		{
			it->second.insert(segl, segr);
			max_coverlength = std::max(max_coverlength, it->second.coverlength());
		}
	}
}

void SegmentCouples::insert(const std::vector<std::pair<unsigned int, unsigned int>>& cp, unsigned int segl, unsigned int segr)
{
	for (auto& p : cp)
	{
		for (unsigned int i = p.first; i <= p.second; i++)
		{
			auto it = storage.find(i);
			if (it == storage.end())
			{
				SegmentCouple sc;
				sc.insert(segl, segr);
				storage[i] = sc;
				max_coverlength = std::max(max_coverlength, storage[i].coverlength());
			}
			else
			{
				it->second.insert(segl, segr);
				max_coverlength = std::max(max_coverlength, it->second.coverlength());
			}
		}
	}
}


std::atomic_bool ErrorRateConventor::m_loglikInited{ false };
std::mutex ErrorRateConventor::m_Mutex;

#ifndef USE_STAT_DATA
std::vector<double> ErrorRateConventor::m_loglik;

ErrorRateConventor::ErrorRateConventor(unsigned int maxlen, unsigned int alphabet_size, unsigned int standard_error)
	: m_maxLen(maxlen), m_alphabetSize(alphabet_size), m_standardError(standard_error)
{
	if (!m_loglikInited.load(std::memory_order_relaxed))
	{
		std::scoped_lock lock(m_Mutex);
		if (!m_loglikInited.load(std::memory_order_acquire))
		{
			m_loglik.resize(maxlen * (maxlen + 1) + 1, -1);
			double logT = log(m_alphabetSize);
			double logF = log(static_cast<double>(m_alphabetSize) / (m_alphabetSize - 1));

			for (unsigned int len = 0; len <= m_maxLen; len++)
			{
				for (unsigned int err = 0; err <= len; err++)
				{
					m_loglik[len * m_maxLen + err] = logC(err, len) + err * logF + (len - err) * logT;
				}
			}

			for (unsigned int len = 0; len <= m_maxLen; len++)
			{
				for (unsigned int err = 1; err <= len; err++)
				{
					double& thisl = m_loglik[len * m_maxLen + err];
					double& lastl = m_loglik[len * m_maxLen + err - 1];
					if (thisl > lastl)
					{
						thisl -= log1p(exp(thisl - lastl));
					}
					else
					{
						thisl = lastl - log1p(exp(lastl - thisl));
					}
				}
			}

			m_loglikInited.store(true, std::memory_order_release);
		}
	}

	stLogP = logP(m_maxLen, m_standardError);
}

double ErrorRateConventor::logP(const unsigned int& len, const unsigned int& err) const
{
	return m_loglik[len * m_maxLen + err];
}

bool ErrorRateConventor::goodMatch(const unsigned int& len, const unsigned int& err) const
{
	if (len * m_maxLen + err >= m_loglik.size())
	{
		std::cout << len << " " << err << std::endl;
		throw std::runtime_error("out of range");
	}
	return m_loglik[len * m_maxLen + err] > stLogP;
}
#else
std::vector<double> ErrorRateConventor::m_minl;

ErrorRateConventor::ErrorRateConventor(unsigned int maxlen, unsigned int alphabet_size, unsigned int standard_error)
{
	if (!m_loglikInited.load(std::memory_order_relaxed))
	{
		std::scoped_lock lock(m_Mutex);
		if (!m_loglikInited.load(std::memory_order_acquire))
		{
			double error_rate = static_cast<double>(standard_error) / static_cast<double>(maxlen);
			unsigned int stdl = static_cast<unsigned int>(floor(error_rate / 0.05));
			unsigned int stdr = static_cast<unsigned int>(ceil(error_rate / 0.05));
			double prop = 0.5;
			if (stdl != stdr)
			{
				prop = stdr - error_rate / 0.05;
			}
			if (stdr >= 12)
			{
				throw std::runtime_error("The error rate is higher than what stat data supports, please undefine value USE_STAT_DATA");
			}
			m_minl.resize(21);
			for (unsigned int i = 0; i < 21; i++)
			{
				m_minl[i] = static_cast<double>(STAT_DATA[stdl][i]) * prop + static_cast<double>(STAT_DATA[stdr][i]) * (1 - prop);
			}
			m_loglikInited.store(true, std::memory_order_release);
		}
	}
}

bool ErrorRateConventor::goodMatch(const unsigned int& len, const unsigned int& err) const
{
	if (err >= m_minl.size())
	{
		std::cout << len << " " << err << std::endl;
		throw std::runtime_error("out of range");
	}
	return len >= m_minl[err];
}

#endif // !USE_STAT_DATA

void ChenAlignerCore::workerThread()
{
	unsigned int task_acquire_each_time = 20;	// the tasks to be transfered to the worker once it has the lock
	// in a new design, it will automatically adjust to pick the tasks that it can finish in around 10s
	unsigned int myTask_beginID;
	bool should_break = false;
	while (true)
	{
		m_Mutex.lock();

		if (m_TasksTaken >= m_TaskReads.size())
		{
			should_break = true;
		}
		else
		{
			std::cout << m_TasksTaken << " Reads Matched" << std::endl;
			myTask_beginID = m_TasksTaken;
			m_TasksTaken += task_acquire_each_time;
		}

		m_Mutex.unlock();

		if (should_break)
			break;

		auto start = std::chrono::high_resolution_clock::now();
		for (unsigned int i = myTask_beginID; i < std::min(myTask_beginID + task_acquire_each_time, static_cast<unsigned int>(m_TaskReads.size())); i++)
		{
			auto ans = alignRead(m_TaskReads[i], m_sub_err, m_all_err, m_err_dif);
			for (auto& a : ans)
			{
				m_Results[i].push_back(a.first);
			}
		}
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

		if (duration != 0)
		{
			task_acquire_each_time = task_acquire_each_time * 10000 / duration;
		}
		if (task_acquire_each_time < 20)
		{
			task_acquire_each_time = 20;
		}
	}
}


ChenAlignerCore::ChenAlignerCore(unsigned int dict_size) : m_CharacterNum(dict_size)
{
	GSANode n0;
	n0.len = 0;
	n0.link = -1;
	m_Nodes.push_back(n0);
	m_IsEraseChannel = false;
}

ChenAlignerCore::ChenAlignerCore(unsigned int dict_size, char erase_marker, unsigned int termination_err_cnt)
	: m_CharacterNum(dict_size), m_EraseMarker(erase_marker), m_TerminationErrorCount(termination_err_cnt)
{
	GSANode n0;
	n0.len = 0;
	n0.link = -1;
	m_Nodes.push_back(n0);
	m_IsEraseChannel = true;
}

bool ChenAlignerCore::addString(const std::string& str, const unsigned int& id)
{
	if (m_HaveBuilt)
		return false;
	if (id == 0)
		throw std::runtime_error("id shouldn't be 0, which is reserved for no matching!");

	m_StringCount++;
	m_LastNode = 0;
	m_StringId = id;
	for (const char& c : str)
	{
		insertTrie(c);
	}
	return true;
}

void ChenAlignerCore::build()
{
	std::cout << "Build starts.\n";
	std::cout << "Step 1: mark watermarks' ID on tire\n";
	buildHostRangeTire();
	std::cout << "Step 2: build suffix automation\n";
	buildSAM();
	std::cout << "Step 3: find parents of nodes\n";
	buildInverseNext();
	std::cout << "Step 4: build tree of inversed links\n";
	buildInverseLink();
	std::cout << "Step 5: mark watermarks' ID on GSA\n";
	buildHostRangeFull();
	if (m_IsEraseChannel)
	{
		std::cout << "(EXTRA) Step 6: calculate erasure link\n";
		buildDeepestLink();
	}
	std::cout << "All steps done.\n";
}

bool ChenAlignerCore::addString(const char* str, const unsigned int& id)
{
	if (m_HaveBuilt)
		return false;
	if (id == 0)
		throw std::runtime_error("id shouldn't be 0, which is reserved for no matching!");

	m_StringCount++;
	m_LastNode = 0;
	m_StringId = id;
	for (unsigned int i = 0; str[i]; i++)
	{
		insertTrie(str[i]);
	}
	return true;
}

void ChenAlignerCore::buildSAM()
{
	m_HaveBuilt = true;
	m_UnclonedLength = m_Nodes.size();

	std::queue<std::pair<unsigned int, char>> awaiting_nodes;	// queue for BFS <parent idx, node idx>
	for (const auto& nd : m_Nodes[0].next)
	{
		awaiting_nodes.push({ 0, nd.first });
	}
	while (!awaiting_nodes.empty())
	{
		unsigned int working_node = insertSAM(awaiting_nodes.front().first, awaiting_nodes.front().second);
		awaiting_nodes.pop();
		for (const auto& nd : m_Nodes[working_node].next)
		{
			awaiting_nodes.push({ working_node, nd.first });
		}
	}
}

void ChenAlignerCore::buildInverseLink()
{
	for (unsigned int i = 1; i < m_Nodes.size(); i++)
	{
		m_Nodes[m_Nodes[i].link].inverse_link.push_back(i);
	}
}

void ChenAlignerCore::buildInverseNext()
{
	m_NodeBefore.resize(m_Nodes.size());
	std::queue<unsigned int> q;
	q.push(0);

	while (!q.empty())
	{
		for (auto& m : m_Nodes[q.front()].next)
		{
			if (m_Nodes[m.second].len == m_Nodes[q.front()].len + 1)
			{
				m_NodeBefore[m.second] = q.front();
				q.push(m.second);
			}
		}
		q.pop();
	}
}

void ChenAlignerCore::calFullHostRange(const unsigned int& node)
{
	SegmentCouple SC;
	if (isNotCloned(node))
	{
		SC.insert(m_HostRangeTrie[node].first, m_HostRangeTrie[node].second);
	}
	for (auto& v : m_Nodes[node].inverse_link)
	{
		calFullHostRange(v);
		for (auto& seg : m_HostRange[v].get_all_segments())
		{
			SC.insert(seg.first, seg.second);
		}
	}
	m_HostRange[node] = SC;
	m_HostRangeCache[node] = SC.get_all_segments();
}

void ChenAlignerCore::buildHostRangeFull()
{
	m_HostRange.resize(m_Nodes.size());
	m_HostRangeCache.resize(m_Nodes.size());
	calFullHostRange(0);
}

void ChenAlignerCore::buildDeepestLink()
{
	m_DeepestLink.resize(m_Nodes.size());
	m_EraseChainLen.resize(m_Nodes.size());

	struct WorkingNode
	{
		bool m_isEraseMarker;
		unsigned int m_Id;
		unsigned int m_BestLinkBefore;
		unsigned int m_ChainCount;

		WorkingNode(bool is_erase, unsigned int id, unsigned int best_before, unsigned int count)
			: m_isEraseMarker(is_erase), m_Id(id), m_BestLinkBefore(best_before), m_ChainCount(count) {}
	};

	std::queue<WorkingNode> node_queue;
	m_DeepestLink[0] = -1;
	for (auto& nd : m_Nodes[0].next)
	{
		node_queue.push(WorkingNode(nd.first == m_EraseMarker, nd.second, 0, nd.first == m_EraseMarker ? 1 : 0));
	}

	while (!node_queue.empty())
	{
		WorkingNode u = std::move(node_queue.front());
		node_queue.pop();

		if (!u.m_isEraseMarker)
		{
			m_DeepestLink[u.m_Id] = m_Nodes[u.m_Id].link;
			m_EraseChainLen[u.m_Id] = 0;
			for (auto& v : m_Nodes[u.m_Id].next)
			{
				if (m_Nodes[v.second].len == m_Nodes[u.m_Id].len + 1)
					node_queue.push(WorkingNode(v.first == m_EraseMarker, v.second, m_Nodes[u.m_Id].link, 0));
			}
		}
		else
		{
			unsigned int best_link = u.m_BestLinkBefore;
			if (m_Nodes[u.m_BestLinkBefore].len < m_Nodes[m_Nodes[u.m_Id].link].len)
			{
				best_link = m_Nodes[u.m_Id].link;
			}
			m_DeepestLink[u.m_Id] = best_link;
			m_EraseChainLen[u.m_Id] = u.m_ChainCount;
			for (auto& v : m_Nodes[u.m_Id].next)
			{
				if (m_Nodes[v.second].len == m_Nodes[u.m_Id].len + 1)
					node_queue.push(WorkingNode(v.first == m_EraseMarker, v.second, best_link, u.m_ChainCount + 1));
			}
		}
	}
}

//void ChenAlignerCore::buildParentMapping()
//{
//	std::stack<unsigned int> s;		// stack for dfs
//	s.push(0);
//
//	while (!s.empty())
//	{
//		unsigned int uidx = s.top();
//		GSANode& u = m_Nodes[uidx];
//		s.pop();
//		for (auto& nxt : u.next)
//		{
//			s.push(nxt.second);
//			GSANode& v = m_Nodes[nxt.second];
//			if (!v.link)
//			{
//				v.parent = nxt.second;
//			}
//			else if (m_Nodes[u.link].next[nxt.first] == v.link)
//			{
//				v.parent = u.parent;
//			}
//			else
//			{
//				v.parent = s.top();
//				// assert (u.len == 1);
//			}
//		}
//	}
//}

void ChenAlignerCore::calHostRange(const unsigned int& node)
{
	GSANode& u = m_Nodes[node];

	if (!u.next.size())
	{
		m_HostRangeTrie[node].first = m_HostMapping.size();
		m_HostRangeTrie[node].second = m_HostMapping.size();
		m_HostMapping.push_back(u.host);
	}
	else
	{
		unsigned int l = m_StringCount, r = 0;
		for (auto& m : u.next)
		{
			calHostRange(m.second);
			l = std::min(l, m_HostRangeTrie[m.second].first);
			r = std::max(r, m_HostRangeTrie[m.second].second);
		}
		m_HostRangeTrie[node].first = l;
		m_HostRangeTrie[node].second = r;
	}
}

std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> ChenAlignerCore::oldFunction_exactCommonSubstring(const std::string& str)
{
	unsigned int t = 0, l = 0, bt = 0, bl = 0;	// current state and matched length, and the best of both
	unsigned int p = 0, bp = 0;					// position on 'str' (end position + 1)

	for (const char& c : str)
	{
		p++;
		while (t && !m_Nodes[t].next[c])
		{
			t = m_Nodes[t].link;
			l = m_Nodes[t].len;
		}
		if (m_Nodes[t].next[c])
		{
			t = m_Nodes[t].next[c];
			l++;
		}
		if (l > bl)
		{
			bt = t;
			bl = l;
			bp = p;
		}
	}

	GSANode& best = m_Nodes[bt];
	return { best.host, best.pos - bl, bp - bl, bl };
}

std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int>> ChenAlignerCore::oldFunction_alignWithSubstitution(const std::string& str, const unsigned int& maxE)
{
	std::priority_queue<MatchedSubstring> longest_ones;
	unsigned int maxSize = 7;		// maximum size of 'longest_ones'

	/// init segment tree
	SegmentTree ST(str.size());

	/// Exact Matching Part
	std::vector<unsigned int> lenFront, lenBack;	// longest suffix of str[0:i] and longest prefix of str[i:-1]
	std::vector<unsigned int> endNode;				// the state after matching the i_th char
	lenFront.resize(str.size());
	lenBack.resize(str.size());
	endNode.resize(str.size());
	unsigned int t = 0, l = 0, p = 0;

	for (const char& c : str)
	{
		if (m_Nodes[t].next[c])
		{
			t = m_Nodes[t].next[c];
			l++;
		}
		else	// a change of beginning will happen
		{
			longest_ones.push({ t, p, l, 0 });
			while (longest_ones.size() > maxSize) longest_ones.pop();
			ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);
			unsigned int lastl = l;
			while (t && !m_Nodes[t].next[c])
			{
				t = m_Nodes[t].link;
				l = m_Nodes[t].len;
			}
			for (; lastl > l; lastl--)
			{
				lenBack[p - lastl] = lastl;
			}
			if (m_Nodes[t].next[c])
			{
				t = m_Nodes[t].next[c];
				l++;
			}
		}
		lenFront[p] = l;
		endNode[p] = t;
		p++;
	}
	longest_ones.push({ t, p, l, 0 });
	while (longest_ones.size() > maxSize) longest_ones.pop();
	ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);
	unsigned int lastl = l;
	for (; lastl > 0; lastl--)
	{
		lenBack[p - lastl] = lastl;
	}

	/// Joining Part
	/* Algorithm
	* Assuming the i_th char of str is a substitution
	* then the upper bound of the length of the substring with only this error
	* is lenFront[i-1] + lenBack[i+1] + 1
	* If this length is shorter than the longest substring containing this char found in less-error match
	* then it can be discarded
	* Otherwise run it on the SAM without link, to see if these two parts are from the same host
	* The longestest string containing one char is maintained by a modified Segment Tree (support segmentary change and elementary query)
	* Error is allowed in lenFront but not in lenBack
	*/
	auto cmp_ep = [](const std::pair<unsigned int, unsigned int>& a, const std::pair<unsigned int, unsigned int>& b) {return a.first < b.first; };
	std::priority_queue<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int>>, decltype(cmp_ep)> error_pos(cmp_ep);

	for (unsigned int err = 1; err <= maxE; err++)
	{
		for (unsigned int i = 1; i < str.size() - 1; i++)
		{
			error_pos.push({ lenFront[i - 1] + lenBack[i + 1] + 1 , i });
		}
		std::vector<unsigned int> lenFrontE = lenFront;
		std::vector<unsigned int> endNodeE = endNode;
		while (!error_pos.empty())
		{
			unsigned int i = error_pos.top().second;
			error_pos.pop();
			unsigned int old_best_len = ST.Query(i, 0, str.size() - 1, 0), best_len = 0;
			//std::cout << i << "," << str[i] << ": " << old_best_len << "->";
			//std::cout << "[" << lenFront[i - 1] << "+" << lenBack[i + 1] << "+1=" << lenFront[i - 1] + lenBack[i + 1] + 1 << "]" << std::endl;
			if (lenFront[i - 1] + lenBack[i + 1] + 1 > old_best_len)
			{
				//std::cout << "\tNode" << endNode[i - 1] << " l=" << m_Nodes[endNode[i - 1]].len << std::endl;
				for (auto& m : m_Nodes[endNode[i - 1]].next)
				{
					if (m.first != str[i])
					{
						unsigned int best_pos, best_node;
						best_len = old_best_len;
						substringExtension(m.second, lenFront[i - 1] + 1, lenBack[i + 1], str, i + 1, best_len, best_node, best_pos);
						if (best_len > old_best_len)
						{
							//std::cout << "\tnew best  [" << best_pos - best_len + 1 << "-" << best_pos << "]->" << best_len << std::endl;
							longest_ones.push({ best_node, best_pos + 1, best_len, err });
							while (longest_ones.size() > maxSize) longest_ones.pop();
							ST.Update(best_pos + 1 - best_len, best_pos, best_len, 0, str.size() - 1, 0);	// not best_pos - 1, because best_pos is updated after that pos matches
							for (; best_pos >= i; best_pos--, best_node = m_NodeBefore[best_node], best_len--)	// >: no continous error; >=: have continous error
							{
								if (best_len > lenFrontE[best_pos])
								{
									lenFrontE[best_pos] = best_len;
									endNode[best_pos] = best_node;
								}
								else
								{
									break;		// the char before will have best length at least this-1
								}
							}
						}
					}
				}
			}
		}
		lenFront = std::move(lenFrontE);
		endNode = std::move(endNodeE);
	}

	std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int>> ans;
	while (!longest_ones.empty())
	{
		const MatchedSubstring& s = longest_ones.top();
		GSANode& n = m_Nodes[s.node_num];
		ans.push_back({ n.host, n.pos - s.matched_len, s.pattern_pos - s.matched_len, s.matched_len, s.error_number });
		longest_ones.pop();
	}
	return ans;
}


std::vector<std::pair<unsigned int, SegmentCouple>> ChenAlignerCore::alignRead(const std::string& str, const unsigned int& sub_err, const unsigned int& all_err, const unsigned int& err_dif) const
{
#ifdef SHOW_OUTPUT_substring_history
	std::ofstream hist_file("len_err_hist/" + channelTypeName + "/" + read2filename(str) + ".len_err.csv");
	std::ofstream sc_hist("len_err_hist/" + channelTypeName + "/" + read2filename(str) + ".seg_cp.csv");
#endif

	if (m_IsEraseChannel)
		return alignReadWithErasure(str, m_EraseMarker, sub_err, all_err, err_dif);

	ErrorRateConventor subJudger(str.size(), m_Alphabet.size(), sub_err);	// to decide whether to add a matched substring to final list. Reject if error density is high.
	SegmentCouples segments;  // stores the accepted substrings

	/// init segment tree
	SegmentTree ST(str.size());  // segment tree to store the longest substring a char at each position on the read belongs to which can be matched with limited hamming distance to any substring of any watermark

	/// Exact Matching Part
	std::vector<unsigned int> lenFront, lenBack;	// longest suffix of str[0:i] and longest prefix of str[i:-1]
	std::vector<unsigned int> endNode;				// the state after matching the i_th char
	std::vector<bool> endNodeWasUpdated;			// useful for ignoring positions that won't change
	lenFront.resize(str.size());
	lenBack.resize(str.size());
	endNode.resize(str.size());
	endNodeWasUpdated.resize(str.size(), true);
	unsigned int t = 0, l = 0, p = 0;

	for (const char& c : str)
	{
		if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
		{
			t = m_Nodes[t].next.at(c);
			l++;
		}
		else	// a change of beginning will happen
		{
			if (subJudger.goodMatch(l, 0))
			{
				segments.insert(m_HostRangeCache[t], p - l, p - 1);
			}

			ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);
			unsigned int lastl = l;
			while (t && m_Nodes[t].next.find(c) == m_Nodes[t].next.end())
			{
				t = m_Nodes[t].link;
				l = m_Nodes[t].len;
			}
			for (; lastl > l; lastl--)
			{
				lenBack[p - lastl] = lastl;
			}
			if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
			{
				t = m_Nodes[t].next.at(c);
				l++;
			}
		}
		lenFront[p] = l;
		endNode[p] = t;
		p++;
	}
#ifdef SHOW_OUTPUT_substring_history
		hist_file << l << "," << "0,";
		print_host_pairs(m_HostRangeCache[t], hist_file);
		hist_file << "\n";
#endif
	if (subJudger.goodMatch(l, 0))
	{
		segments.insert(m_HostRangeCache[t], p - l, p - 1);
	}
	ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);
	unsigned int lastl = l;
	for (; lastl > 0; lastl--)
	{
		lenBack[p - lastl] = lastl;
	}

	/// Joining Part
	auto cmp_ep = [](const std::pair<unsigned int, unsigned int>& a, const std::pair<unsigned int, unsigned int>& b) {return a.first < b.first; };
	std::priority_queue<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int>>, decltype(cmp_ep)> error_pos(cmp_ep);

	unsigned int max_loop_err = std::min(m_TerminationErrorCount, sub_err);
	for (unsigned int err = 1; err <= max_loop_err && segments.maxcoverlength() + all_err < str.size(); err++)  // breadth first searching, gradually adding errors
	{
#ifdef SHOW_OUTPUT_substring_history
		sc_hist << segments.storage.size() << "," << segments.max_coverlength << "\n";
#endif
		ST.updateTagReset();
		for (unsigned int i = 1; i < str.size() - 1; i++)  // sort the substrings and start matching from the most promising one
		{
			if (endNodeWasUpdated[i])
			{
				error_pos.push({ lenFront[i - 1] + lenBack[i + 1] + 1 , i });
			}
		}
		if (error_pos.empty())
		{
			break;
		}
		std::vector<unsigned int> lenFrontE = lenFront;
		std::vector<unsigned int> endNodeE = endNode;
		while (!error_pos.empty())
		{
			unsigned int i = error_pos.top().second;
			error_pos.pop();
			unsigned int old_best_len = ST.Query(i, 0, str.size() - 1, 0), best_len = 0;
			if (lenFront[i - 1] + lenBack[i + 1] + 1 > old_best_len)
				// cancel searching if the upper bound of this searching is beaten by another substring
			{
				for (auto& m : m_Nodes[endNode[i - 1]].next)  // to search the grandchildren of the node before the mismatch, this loop first list all its children
				{
					if (m.first != str[i])  // we assumed this position is an error, and there's no need to assume a match is an error (waste of time)
					{
						unsigned int best_pos = i;
						unsigned int best_node = m.second;
						best_len = std::max(old_best_len, lenFront[i - 1] + 1);		// the second term is usually smaller except at the end of one matched substring
						auto len_before = substringExtension(m.second, lenFront[i - 1] + 1, lenBack[i + 1], str, i + 1, best_len, best_node, best_pos);
						if (best_len > old_best_len)  // reject if extending from here can't reach a better result
						{
#ifdef SHOW_OUTPUT_substring_history
							hist_file << best_len << "," << err << ",";
							print_host_pairs(m_HostRangeCache[best_node], hist_file);
							hist_file << "\n";
#endif
							if (subJudger.goodMatch(best_len, err))  // insert good substring match to final list
							{
								segments.insert(m_HostRangeCache[best_node], best_pos + 1 - best_len, best_pos);
							}
							// update BestMatchBelong
							ST.Update(best_pos + 1 - best_len, best_pos, best_len, 0, str.size() - 1, 0);	// not best_pos - 1, because best_pos is updated after that pos matches
							// update BestMatchFront
							for (; best_pos >= i; best_pos--, best_len--)	// >: no continous error; >=: have continous error
							{
								if (best_len > lenFrontE[best_pos])
								{
									lenFrontE[best_pos] = best_len;
									endNode[best_pos] = best_node;
								}
								else
								{
									break;		// the char before will have best length at least this-1
								}

								best_node = m_NodeBefore[best_node];
								if (len_before.size() <= 1)
								{
									continue;
								}
								len_before.pop_back();
								while (best_node && m_Nodes[best_node].len > len_before.back())
								{
									best_node = m_Nodes[best_node].link;
								}
								if (m_Nodes[best_node].len < len_before.back())
								{
									throw std::runtime_error("IMPOSSIBLE?! m_Nodes[best_node].len < len_before.back()");
								}
							}
							// note that BestMatchEnd is always error free and don't need update
						}
					}
				}
			}
		}
		lenFront = std::move(lenFrontE);  // O(0) move
		for (unsigned int i = 0; i < endNode.size(); i++)
		{
			endNodeWasUpdated[i] = endNode[i] == endNodeE[i] ? false : true;
		}
		endNode = std::move(endNodeE);
	}

#ifdef SHOW_OUTPUT_substring_history
	hist_file.close();
	sc_hist.close();
#endif
	return segments.show_matches(err_dif, m_HostMapping);
}

inline void ChenAlignerCore::print_host_pairs(const std::vector<std::pair<unsigned int, unsigned int>> data, std::ostream& os) const
{
	for (auto& p : data)
	{
		os << p.first << "," << p.second << ",";
	}
}

std::vector<std::pair<unsigned int, SegmentCouple>> ChenAlignerCore::alignReadWithErasure(const std::string& str, const char& erased_char, const unsigned int& nature_sub_err, const unsigned int& nature_other_err, const unsigned int& err_dif) const
{
#ifdef SHOW_OUTPUT_substring_history
	std::ofstream hist_file("len_err_hist/" + channelTypeName + "/" + read2filename(str) + ".len_err.csv");
	std::ofstream sc_hist("len_err_hist/" + channelTypeName + "/" + read2filename(str) + ".seg_cp.csv");
#endif

	ErrorRateConventor subJudger(str.size(), m_Alphabet.size(), nature_sub_err);	// to decide whether to add a matched substring to final list. Reject if error density is high.
	SegmentCouples segments;  // stores the accepted substrings

	/// init segment tree
	SegmentTree ST(str.size());  // segment tree to store the longest substring a char at each position on the read belongs to which can be matched with limited hamming distance to any substring of any watermark


	/// Exact Matching Part
	std::vector<unsigned int> lenFront, lenBack;	// longest suffix of str[0:i] and longest prefix of str[i:-1]
	std::vector<unsigned int> endNode;				// the state after matching the i_th char
	std::vector<bool> endNodeWasUpdated;			// useful for ignoring positions that won't change
	lenFront.resize(str.size());
	lenBack.resize(str.size());
	endNode.resize(str.size());
	endNodeWasUpdated.resize(str.size(), true);
	unsigned int t = 0, l = 0, p = 0;

	if constexpr (INIT_TYPE == 0)
	{
		for (const char& c : str)
		{
			if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
			{
				t = m_Nodes[t].next.at(c);
				l++;
			}
			else	// a change of beginning will happen
			{
#ifdef SHOW_OUTPUT_substring_history
				hist_file << l << "," << "0,";
				print_host_pairs(m_HostRangeCache[t], hist_file);
				hist_file << "\n";
#endif
				if (subJudger.goodMatch(l, 0))
				{
					segments.insert(m_HostRangeCache[t], p - l, p - 1);
					if constexpr (SHOW_OUTPUT_print_segment)
					{
						print_segment(str, m_HostRangeCache[t], p - l, p - 1);
					}
				}

				ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);
				unsigned int lastl = l;
				while (t && m_Nodes[t].next.find(c) == m_Nodes[t].next.end())
				{
					t = m_Nodes[t].link;
					l = m_Nodes[t].len;
				}
				for (; lastl > l; lastl--)
				{
					lenBack[p - lastl] = lastl;
				}
				if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
				{
					t = m_Nodes[t].next.at(c);
					l++;
				}
			}
			lenFront[p] = l;
			endNode[p] = t;
			p++;
		}
#ifdef SHOW_OUTPUT_substring_history
		hist_file << l << "," << "0,";
		print_host_pairs(m_HostRangeCache[t], hist_file);
		hist_file << "\n";
#endif
		if (subJudger.goodMatch(l, 0))
		{
			segments.insert(m_HostRangeCache[t], p - l, p - 1);
			if constexpr (SHOW_OUTPUT_print_segment)
			{
				print_segment(str, m_HostRangeCache[t], p - l, p - 1);
			}
		}
		ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);
		unsigned int lastl = l;
		for (; lastl > 0; lastl--)
		{
			lenBack[p - lastl] = lastl;
		}
	}
	else if constexpr (INIT_TYPE == 1)
	{
		/* This part is discarded, because it can't accurately estimate the lenBack.
		A more accurate but (perhaps) slower algorithm is used for replacement.*/

		// best back use maximum match including any erasure
		// best front must end with a non-erasure character

		// don't consider erasure and update lenFront
		for (const char& c : str)
		{
			if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
			{
				t = m_Nodes[t].next.at(c);
				l++;
			}
			else	// a change of beginning will happen
			{
				while (t && m_Nodes[t].next.find(c) == m_Nodes[t].next.end())
				{
					t = m_Nodes[t].link;
					l = m_Nodes[t].len;
				}
				if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
				{
					t = m_Nodes[t].next.at(c);
					l++;
				}
			}
			lenFront[p] = l;
			endNode[p] = t;
			p++;
		}

		// consider erasure and update all, inclusing lenFront
		t = 0; l = 0; p = 0;
		unsigned int best_link_depth = 0, best_link_source, best_link_p;
		for (p = 0; p < str.size(); p++)
		{
			char c = str[p];
			if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
			{
				t = m_Nodes[t].next.at(c);
				l++;

				lenFront[p] = l;	// changed place of update
				endNode[p] = t;

				best_link_depth = m_Nodes[m_Nodes[t].link].len;
				best_link_source = t;
				best_link_p = p;
			}
			else if (m_Nodes[t].next.find(m_EraseMarker) != m_Nodes[t].next.end())
			{
				t = m_Nodes[t].next.at(m_EraseMarker);
				l++;

				if (m_Nodes[m_Nodes[t].link].len >= best_link_depth)
				{
					best_link_depth = m_Nodes[m_Nodes[t].link].len;
					best_link_source = t;
					best_link_p = p;
				}
			}
			else	// fail to search
			{
#ifdef SHOW_OUTPUT_substring_history
				hist_file << l << "," << "0,";
				print_host_pairs(m_HostRangeCache[t], hist_file);
				hist_file << "\n";
#endif
				if (subJudger.goodMatch(l, 0))
				{
					segments.insert(m_HostRangeCache[t], p - l, p - 1);
					if constexpr (SHOW_OUTPUT_print_segment)
					{
						print_segment(str, m_HostRangeCache[t], p - l, p - 1);
					}
				}

				ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);
				unsigned int lastl = l;
				for (; lastl > l; lastl--)
				{
					lenBack[p - lastl] = lastl;
				}

				if (best_link_depth > 0)
				{
					best_link_depth = 0;
					p = best_link_p + 1;
					c = str[p];
					t = m_Nodes[best_link_source].link;
					l = m_Nodes[t].len;
				}
				unsigned int original_depth = m_EraseChainLen[t];
				while (t && m_Nodes[t].next.find(c) == m_Nodes[t].next.end()
					&& m_Nodes[t].next.find(m_EraseMarker) == m_Nodes[t].next.end())
				{
					t = m_DeepestLink[t];
					l = m_Nodes[t].len;
					p = p + m_EraseChainLen[t] - original_depth;
					c = str[p];
					original_depth = m_EraseChainLen[t];
				}
			}
		}
#ifdef SHOW_OUTPUT_substring_history
		hist_file << l << "," << "0,";
		print_host_pairs(m_HostRangeCache[t], hist_file);
		hist_file << "\n";
#endif
		if (subJudger.goodMatch(l, 0))
		{
			segments.insert(m_HostRangeCache[t], p - l, p - 1);
			if constexpr (SHOW_OUTPUT_print_segment)
			{
				print_segment(str, m_HostRangeCache[t], p - l, p - 1);
			}
		}
		ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);
		unsigned int lastl = l;
		for (; lastl > 0; lastl--)
		{
			lenBack[p - lastl] = lastl;
		}
	}
	else if constexpr (INIT_TYPE == 2)
	{
		std::vector<std::set<unsigned int>> visited_nodes;
		visited_nodes.resize(str.size());
		for (unsigned int p0 = 0; p0 < str.size(); p0++)
		{
			t = 0; l = 0;
			for (p = p0; p < str.size(); p++)
			{
				const char& c = str[p];
				if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
				{
					t = m_Nodes[t].next.at(c);
					l++;
				}
				else if (m_Nodes[t].next.find(m_EraseMarker) != m_Nodes[t].next.end())
				{
					t = m_Nodes[t].next.at(m_EraseMarker);
					l++;
				}
				else	// a change of beginning will happen
				{
					break;
				}

				if (visited_nodes[p].count(t))
				{
					break;
				}
				visited_nodes[p].insert(t);
				if (l > lenFront[p])
				{
					lenFront[p] = l;
					endNode[p] = t;
				}
			}
		}

		for (int i = str.size() - 1; i >= 0; i--)
		{
			for (int j = lenFront[i]; j > 0; j--)
			{
				if (lenBack[i - j + 1] >= j)
				{
					break;
				}
				else
				{
					lenBack[i - j + 1] = j;
				}
			}
		}
	}
	else if constexpr (INIT_TYPE == 3)
	{
		// use a broader len back
		for (const char& c : str)
		{
			if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
			{
				t = m_Nodes[t].next.at(c);
				l++;
			}
			else	// a change of beginning will happen
			{
#ifdef SHOW_OUTPUT_substring_history
				hist_file << l << "," << "0,";
				print_host_pairs(m_HostRangeCache[t], hist_file);
				hist_file << "\n";
#endif
				if (subJudger.goodMatch(l, 0))
				{
					segments.insert(m_HostRangeCache[t], p - l, p - 1);
					if constexpr (SHOW_OUTPUT_print_segment)
					{
						print_segment(str, m_HostRangeCache[t], p - l, p - 1);
					}
				}

				ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);
				unsigned int lastl = l;
				while (t && m_Nodes[t].next.find(c) == m_Nodes[t].next.end())
				{
					t = m_Nodes[t].link;
					l = m_Nodes[t].len;
				}
				if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
				{
					t = m_Nodes[t].next.at(c);
					l++;
				}
			}
			lenFront[p] = l;
			endNode[p] = t;
			p++;
		}
#ifdef SHOW_OUTPUT_substring_history
		hist_file << l << "," << "0,";
		print_host_pairs(m_HostRangeCache[t], hist_file);
		hist_file << "\n";
#endif
		if (subJudger.goodMatch(l, 0))
		{
			segments.insert(m_HostRangeCache[t], p - l, p - 1);
			if constexpr (SHOW_OUTPUT_print_segment)
			{
				print_segment(str, m_HostRangeCache[t], p - l, p - 1);
			}
		}
		ST.Update(p - l, p - 1, l, 0, str.size() - 1, 0);

		for (unsigned int i = 0; i < lenBack.size(); i++)
		{
			lenBack[i] = str.size() - i;
		}
	}

	/// Joining Part
	auto cmp_ep = [](const std::pair<unsigned int, unsigned int>& a, const std::pair<unsigned int, unsigned int>& b) {return a.first < b.first; };
	std::priority_queue<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int>>, decltype(cmp_ep)> error_pos(cmp_ep);

	unsigned int max_loop_err = std::min(m_TerminationErrorCount, nature_sub_err);
	for (unsigned int err = 1; err <= max_loop_err && segments.maxcoverlength() + nature_other_err < str.size(); err++)  // breadth first searching, gradually adding errors
	{
#ifdef SHOW_OUTPUT_substring_history
		sc_hist << segments.storage.size() << "," << segments.max_coverlength << "\n";
#endif
		ST.updateTagReset();
		for (unsigned int i = 1; i < str.size() - 1; i++)  // sort the substrings and start matching from the most promising one
		{
			if (endNodeWasUpdated[i])
			{
				error_pos.push({ lenFront[i - 1] + lenBack[i + 1] + 1 , i });
			}
		}
		if (error_pos.empty())
		{
			break;
		}
		std::vector<unsigned int> lenFrontE = lenFront;
		std::vector<unsigned int> endNodeE = endNode;
		while (!error_pos.empty())
		{
			unsigned int i = error_pos.top().second;
			error_pos.pop();
			unsigned int old_best_len = ST.Query(i, 0, str.size() - 1, 0), best_len = 0;
			if (lenFront[i - 1] + lenBack[i + 1] + 1 > old_best_len)
				// cancel searching if the upper bound of this searching is beaten by another substring
			{
				for (auto& m : m_Nodes[endNode[i - 1]].next)  // to search the grandchildren of the node before the mismatch, this loop first list all its children
				{
					if (m.first != str[i])  // we assumed this position is an error, and there's no need to assume a match is an error (waste of time)
					{
						unsigned int best_pos = i;
						unsigned int best_node = m.second;
						best_len = std::max(old_best_len, lenFront[i - 1] + 1);		// the second term is usually smaller except at the end of one matched substring
						auto len_before = substringExtensionWithErasure(m.second, lenFront[i - 1] + 1, lenBack[i + 1], str, i + 1, best_len, best_node, best_pos, erased_char);
						if (best_len > old_best_len)  // reject if extending from here can't reach a better result
						{
							if constexpr (SHOW_OUTPUT_enlengthen_procedure)
							{
								for (auto& ll : len_before)
									std::cout << ll << " ";
								std::cout << std::endl;
							}
#ifdef SHOW_OUTPUT_substring_history
							hist_file << best_len << "," << err << ",";
							print_host_pairs(m_HostRangeCache[best_node], hist_file);
							hist_file << "\n";
#endif
							if (subJudger.goodMatch(best_len, err))  // insert good substring match to final list
							{
								segments.insert(m_HostRangeCache[best_node], best_pos + 1 - best_len, best_pos);
								if constexpr (SHOW_OUTPUT_print_segment)
								{
									print_segment(str, m_HostRangeCache[best_node], best_pos + 1 - best_len, best_pos);
								}
							}
							// update BestMatchBelong
							ST.Update(best_pos + 1 - best_len, best_pos, best_len, 0, str.size() - 1, 0);	// not best_pos - 1, because best_pos is updated after that pos matches
							// update BestMatchFront
							std::vector<unsigned int> my_len_hist;
							my_len_hist.push_back(m_Nodes[best_node].len);
							for (; best_pos >= i; best_pos--, best_len--)	// >: no continous error; >=: have continous error
							{
								if (best_len > lenFrontE[best_pos])
								{
									lenFrontE[best_pos] = best_len;
									endNode[best_pos] = best_node;
								}
								else
								{
									break;		// the char before will have best length at least this-1
								}

								my_len_hist.push_back(0);
								best_node = m_NodeBefore[best_node];
								my_len_hist.push_back(m_Nodes[best_node].len);
								if (len_before.size() <= 1)
								{
									continue;
								}
								len_before.pop_back();
								while (best_node && m_Nodes[best_node].len > len_before.back())
								{
									my_len_hist.push_back(m_Nodes[best_node].len);
									best_node = m_Nodes[best_node].link;
								}
								if (m_Nodes[best_node].len < len_before.back())
								{
									for (auto& lh : my_len_hist)
										std::cout << lh << " ";
									std::cout << std::endl;
									std::cout << len_before.size() << " " << len_before.back() << " " << best_node << " " << m_Nodes[best_node].len << std::endl;
									throw std::runtime_error("IMPOSSIBLE?! m_Nodes[best_node].len < len_before.back()");
								}
							}
							// note that BestMatchEnd is always error free and don't need update
						}
					}
				}
			}
		}
		lenFront = std::move(lenFrontE);  // O(0) move

		for (unsigned int i = 0; i < endNode.size(); i++)
		{
			endNodeWasUpdated[i] = endNode[i] == endNodeE[i] ? false : true;
		}

		endNode = std::move(endNodeE);
	}
	
#ifdef SHOW_OUTPUT_substring_history
	hist_file.close();
	sc_hist.close();
#endif

	return segments.show_matches(err_dif, m_HostMapping);
}

void ChenAlignerCore::initThreadPool(unsigned int thread_num)
{
	for (unsigned int i = 0; i < thread_num; i++)
	{
		m_Matchers.emplace_back(std::thread(&ChenAlignerCore::workerThread, this));
	}
}


void ChenAlignerCore::parallel_Init(const unsigned int& sub_err, const unsigned int& all_err, const unsigned int& err_dif)
{
	m_sub_err = sub_err;
	m_all_err = all_err;
	m_err_dif = err_dif;
}

void ChenAlignerCore::parallel_addRead(const std::string& str)
{
	m_TaskReads.push_back(str);
}

void ChenAlignerCore::parallel_addRead(const std::string& str, const unsigned int& sub_err, const unsigned int& all_err, const unsigned int& err_dif)
{
	parallel_addRead(str);
	parallel_Init(sub_err, all_err, err_dif);
}

void ChenAlignerCore::parallel_terminateReadTransfer(const unsigned int& thread_num)
{
	m_Results.resize(m_TaskReads.size());
	m_TasksTaken = 0;

	initThreadPool(thread_num);

	for (auto& t : m_Matchers)
	{
		t.join();
	}
}

void ChenAlignerCore::parallel_getResults(const std::string& outfile_name, OUTFILE_FORMAT format)
{
	std::ofstream file(outfile_name.c_str());
	switch (format)
	{
	case OUTFILE_FORMAT::BEST:
		for (auto& res : m_Results)
		{
			if (res.empty())
				file << 0 << "\n";
			else
				file << res[0] << "\n";
		}
		break;
	case OUTFILE_FORMAT::ALL:
		for (auto& res : m_Results)
		{
			file << res.size() << ",";
			for (auto& i : res)
			{
				file << i << ",";
			}
			file << "\n";
		}
		break;
	default:
		break;
	}
}

void ChenAlignerCore::parallel_getResults(std::vector<unsigned int>& best_result)
{
	for (auto& res : m_Results)
	{
		if (res.empty())
			best_result.push_back(0);
		else
			best_result.push_back(res[0]);
	}
}

void ChenAlignerCore::parallel_getResults(std::vector<std::vector<unsigned int>>& all_results)
{
	all_results = m_Results;
	for (auto& res : all_results)
	{
		if (res.empty())
		{
			res.push_back(0);
		}
	}
}

void GSANode::print(std::ostream& os) const
{
	os.write(reinterpret_cast<const char*>(&len), sizeof(len));
	os.write(reinterpret_cast<const char*>(&link), sizeof(link));

	unsigned int next_size = next.size();
	os.write(reinterpret_cast<const char*>(&next_size), sizeof(next_size));
	for (const auto& child : next)
	{
		os.write(&child.first, sizeof(child.first));
		os.write(reinterpret_cast<const char*>(&child.second), sizeof(child.second));
	}

	os.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
	os.write(reinterpret_cast<const char*>(&host), sizeof(host));

	unsigned int inverseLink_size = inverse_link.size();
	os.write(reinterpret_cast<const char*>(&inverseLink_size), sizeof(inverseLink_size));
	for (unsigned int elem : inverse_link)
	{
		os.write(reinterpret_cast<const char*>(&elem), sizeof(elem));
	}
}

void GSANode::read(std::istream& is)
{
	is.read(reinterpret_cast<char*>(&len), sizeof(len));
	is.read(reinterpret_cast<char*>(&link), sizeof(link));

	unsigned int next_size;
	is.read(reinterpret_cast<char*>(&next_size), sizeof(next_size));
	next.clear();
	for (unsigned int i = 0; i < next_size; i++)
	{
		char c;
		unsigned int nodeID;
		is.read(&c, sizeof(c));
		is.read(reinterpret_cast<char*>(&nodeID), sizeof(nodeID));
		next[c] = nodeID;
	}

	is.read(reinterpret_cast<char*>(&pos), sizeof(pos));
	is.read(reinterpret_cast<char*>(&host), sizeof(host));

	unsigned int inverseLink_size;
	is.read(reinterpret_cast<char*>(&inverseLink_size), sizeof(inverseLink_size));
	inverse_link.clear();
	inverse_link.resize(inverseLink_size);
	for (unsigned int& elem : inverse_link)
	{
		is.read(reinterpret_cast<char*>(&elem), sizeof(elem));
	}
}

void ChenAlignerCore::cache_storeStructure(std::ostream& os) const
{
	os.write(reinterpret_cast<const char*>(&m_CharacterNum), sizeof(m_CharacterNum));
	os.write(reinterpret_cast<const char*>(&m_StringCount), sizeof(m_StringCount));

	unsigned int alphabet_size = m_Alphabet.size();
	os.write(reinterpret_cast<const char*>(&alphabet_size), sizeof(alphabet_size));
	for (char c : m_Alphabet)
	{
		os.write(&c, sizeof(c));
	}

	unsigned int nodes_size = m_Nodes.size();
	os.write(reinterpret_cast<const char*>(&nodes_size), sizeof(nodes_size));
	for (const auto& node : m_Nodes)
	{
		node.print(os);
	}

	unsigned int hostRange_size = m_HostRangeCache.size();
	os.write(reinterpret_cast<const char*>(&hostRange_size), sizeof(hostRange_size));
	for (const auto& pairs : m_HostRangeCache)
	{
		unsigned int pairs_size = pairs.size();
		os.write(reinterpret_cast<const char*>(&pairs_size), sizeof(pairs_size));
		for (const auto& pair : pairs)
		{
			os.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
			os.write(reinterpret_cast<const char*>(&pair.second), sizeof(pair.second));
		}
	}

	unsigned int hostMapping_size = m_HostMapping.size();
	os.write(reinterpret_cast<const char*>(&hostMapping_size), sizeof(hostMapping_size));
	for (auto& hostId : m_HostMapping)
	{
		os.write(reinterpret_cast<const char*>(&hostId), sizeof(hostId));
	}

	unsigned int nodeBefore_size = m_NodeBefore.size();
	os.write(reinterpret_cast<const char*>(&nodeBefore_size), sizeof(nodeBefore_size));
	for (auto& nodeId : m_NodeBefore)
	{
		os.write(reinterpret_cast<const char*>(&nodeId), sizeof(nodeId));
	}

	unsigned int deepestLink_size = m_DeepestLink.size();
	os.write(reinterpret_cast<const char*>(&deepestLink_size), sizeof(deepestLink_size));
	for (auto& link : m_DeepestLink)
	{
		os.write(reinterpret_cast<const char*>(&link), sizeof(link));
	}

	unsigned int chainCnt_size = m_EraseChainLen.size();
	os.write(reinterpret_cast<const char*>(&chainCnt_size), sizeof(chainCnt_size));
	for (auto& cnt : m_EraseChainLen)
	{
		os.write(reinterpret_cast<const char*>(&cnt), sizeof(cnt));
	}
}

void ChenAlignerCore::cache_acquireStructure(std::istream& is)
{
	is.read(reinterpret_cast<char*>(&m_CharacterNum), sizeof(m_CharacterNum));
	is.read(reinterpret_cast<char*>(&m_StringCount), sizeof(m_StringCount));

	// Deserialize the alphabet set
	unsigned int alphabetSize;
	is.read(reinterpret_cast<char*>(&alphabetSize), sizeof(alphabetSize));
	m_Alphabet.clear();
	for (unsigned int i = 0; i < alphabetSize; i++)
	{
		char c;
		is.read(&c, sizeof(c));
		m_Alphabet.insert(c);
	}

	// Deserialize the nodes vector
	unsigned int nodesSize;
	is.read(reinterpret_cast<char*>(&nodesSize), sizeof(nodesSize));
	m_Nodes.clear();
	m_Nodes.resize(nodesSize);
	for (auto& node : m_Nodes)
	{
		node.read(is);
	}

	// Deserialize the host range vector
	unsigned int hostRangeSize;
	is.read(reinterpret_cast<char*>(&hostRangeSize), sizeof(hostRangeSize));
	m_HostRangeCache.clear();
	m_HostRangeCache.resize(hostRangeSize);
	for (auto& pairs : m_HostRangeCache)
	{
		unsigned int pairsSize;
		is.read(reinterpret_cast<char*>(&pairsSize), sizeof(pairsSize));
		pairs.clear();
		pairs.resize(pairsSize);
		for (auto& pair : pairs)
		{
			is.read(reinterpret_cast<char*>(&pair.first), sizeof(pair.first));
			is.read(reinterpret_cast<char*>(&pair.second), sizeof(pair.second));
		}
	}

	// Deserialize the host mapping vector
	unsigned int hostMappingSize;
	is.read(reinterpret_cast<char*>(&hostMappingSize), sizeof(hostMappingSize));
	m_HostMapping.clear();
	m_HostMapping.resize(hostMappingSize);
	for (auto& hostId : m_HostMapping)
	{
		is.read(reinterpret_cast<char*>(&hostId), sizeof(hostId));
	}

	// Deserialize the node before vector
	unsigned int nodeBeforeSize;
	is.read(reinterpret_cast<char*>(&nodeBeforeSize), sizeof(nodeBeforeSize));
	m_NodeBefore.clear();
	m_NodeBefore.resize(nodeBeforeSize);
	for (auto& nodeId : m_NodeBefore)
	{
		is.read(reinterpret_cast<char*>(&nodeId), sizeof(nodeId));
	}

	if (!m_IsEraseChannel)
		return;

	// Deserialize the deepest link vector
	unsigned int deepestLinkSize;
	is.read(reinterpret_cast<char*>(&deepestLinkSize), sizeof(deepestLinkSize));
	m_DeepestLink.clear();
	m_DeepestLink.resize(deepestLinkSize);
	for (auto& link : m_DeepestLink)
	{
		is.read(reinterpret_cast<char*>(&link), sizeof(link));
	}

	// Deserialize the erase chain length vector
	unsigned int chainLenSize;
	is.read(reinterpret_cast<char*>(&chainLenSize), sizeof(chainLenSize));
	m_EraseChainLen.clear();
	m_EraseChainLen.resize(deepestLinkSize);
	for (auto& cnt : m_EraseChainLen)
	{
		is.read(reinterpret_cast<char*>(&cnt), sizeof(cnt));
	}
}

std::pair<unsigned int, unsigned int> ChenAlignerCore::oldFunction_exactAlignment(const std::string& str)
{
	unsigned int t = 0;	// current node / state
	for (const char& c : str)
	{
		if (m_Nodes[t].next[c])
		{
			t = m_Nodes[t].next[c];
		}
		else
		{
			return { 0, 0 };
		}
	}
	return { m_Nodes[t].host, m_Nodes[t].pos - str.size() };
}

void ChenAlignerCore::buildHostRangeTire()
{
	m_HostRangeTrie.resize(m_Nodes.size());
	calHostRange(0);
}

void SegmentTree::Update(const unsigned int& left, const unsigned int& right, const unsigned int& value, const unsigned int& segleft, const unsigned int& segright, const unsigned int& node)
{
	m_updateHappened = true;

	if (value <= m_STMin[node])
	{
		return;
	}

	if (left <= segleft && segright <= right)
	{
		m_STMin[node] = value;
		return;
	}

	unsigned int segmid = segleft + ((segright - segleft) >> 1);
	if (left <= segmid)
	{
		Update(left, right, value, segleft, segmid, stlchild(node));
	}
	if (right > segmid)
	{
		Update(left, right, value, segmid + 1, segright, strchild(node));
	}
	m_STMin[node] = std::min(m_STMin[stlchild(node)], m_STMin[strchild(node)]);		// must have both children, otherwise left <= segleft = segright <= right
}

unsigned int SegmentTree::Query(const unsigned int& x, const unsigned int& segleft, const unsigned int& segright, const unsigned int& node)
{
	// brute force for debugging purpose
	//return m_STMin[x];

	if (segleft == segright)
	{
		return m_STMin[node];
	}

	unsigned int segmid = segleft + ((segright - segleft) >> 1);
	if (x <= segmid)
	{
		if (m_STMin[stlchild(node)] < m_STMin[node])
		{
			return m_STMin[node];
		}
		else
		{
			return Query(x, segleft, segmid, stlchild(node));
		}
	}
	else
	{
		if (m_STMin[strchild(node)] < m_STMin[node])
		{
			return m_STMin[node];
		}
		else
		{
			return Query(x, segmid + 1, segright, strchild(node));
		}
	}
}

inline void ChenAlignerCore::insertTrie(const char& c)
{
	m_Alphabet.insert(c);
	if (m_Nodes[m_LastNode].next[c])
	{
		m_LastNode = m_Nodes[m_LastNode].next[c];
	}
	else
	{
		m_Nodes[m_LastNode].next[c] = m_Nodes.size();
		GSANode cur;
		cur.len = 0;	// use len in building for "visited" function, so set to 0 at first
		cur.link = 0;
		cur.host = m_StringId;
		cur.pos = m_Nodes[m_LastNode].pos + 1;
		m_LastNode = m_Nodes.size();
		m_Nodes.push_back(cur);
	}
}

inline unsigned int ChenAlignerCore::insertSAM(const unsigned int& last, const char& c)
{
	unsigned int cur = m_Nodes[last].next[c];	// currently working node
	GSANode& cur_node = m_Nodes[cur];

	// if this node had been updated, directly return (move to it, making no adjustments)
	if (cur_node.len)
	{
		return cur;
	}
	cur_node.len = m_Nodes[last].len + 1;

	int p;	// suffixes of the last node, and see if they're still the suffix with an additional 'c'
	for (p = m_Nodes[last].link; p != -1; p = m_Nodes[p].link)
	{
		if (!m_Nodes[p].next[c])
		{
			m_Nodes[p].next[c] = cur;
		}
		else
		{
			break;
		}
	}
	if (p == -1)
	{
		cur_node.link = 0;
		return cur;
	}

	// no gap between p and q means that with addition 'c' there still form a suffix, it can be directly used
	int q = m_Nodes[p].next[c];
	if (m_Nodes[p].len + 1 == m_Nodes[q].len)
	{
		cur_node.link = q;
		return cur;
	}

	// otherwise there needs a cloned node of q, representing this connection
	// only copy the nodes that have been visited
	// then replace q
	unsigned int clone = m_Nodes.size();
	GSANode clone_node;
	clone_node.len = m_Nodes[p].len + 1;
	clone_node.link = m_Nodes[q].link;
	clone_node.pos = m_Nodes[q].pos;
	clone_node.host = m_Nodes[q].host;
	for (const auto& i : m_Nodes[q].next)
	{
		if (m_Nodes[i.second].len)
		{
			clone_node.next[i.first] = i.second;
		}
	}
	for (; p != -1 && m_Nodes[p].next[c] == q; p = m_Nodes[p].link)
	{
		m_Nodes[p].next[c] = clone;
	}
	m_Nodes.push_back(clone_node);	/// !!! After this step, cur_node should be discarded, because there might be resize of m_Nodes !!!

	m_Nodes[cur].link = clone;
	m_Nodes[q].link = clone;
	return cur;
}

inline std::vector<unsigned int> ChenAlignerCore::substringExtension(const unsigned int& start_node, const unsigned int& start_len, const unsigned int& max_successor,
	const std::string& str, const unsigned int start_pos, unsigned int& best_len, unsigned int& best_node, unsigned int& best_pos) const
{
	unsigned int t = start_node, l = start_len, ms = max_successor, p = start_pos;
	std::vector<unsigned int> len_history;

	for (; p < str.size() && ms >= 0; p++, ms--)
	{
		const char& c = str[p];
		while (t && m_Nodes[t].next.find(c) == m_Nodes[t].next.end() && l + ms > best_len)
		{
			t = m_Nodes[t].link;
			l = m_Nodes[t].len;
			len_history.clear();
		}

		if (l + ms <= best_len)
		{
			break;
		}

		if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
		{
			t = m_Nodes[t].next.at(c);
			l++;
			len_history.push_back(m_Nodes[t].len);	// l <= m_Nodes[t].len. The matched string is a suffix of the longest string to this node.

			// putting this if in if makes sure that the end char is a match
			if (l > best_len)
			{
				best_node = t;
				best_len = l;
				best_pos = p;
			}
		}
	}
	return len_history;		// Technically can autamatically apply std::move()
}

inline std::vector<unsigned int> ChenAlignerCore::substringExtensionWithErasure(const unsigned int& start_node, const unsigned int& start_len, const unsigned int& max_successor, const std::string& str, const unsigned int start_pos, unsigned int& best_len, unsigned int& best_node, unsigned int& best_pos, const char& ignored_char) const
{
	if constexpr (!SHOW_OUTPUT_enlengthenSearch_ignoreErasure)
	{
		unsigned int t = start_node, l = start_len, ms = max_successor, p = start_pos;
		unsigned int best_link_depth = 0, best_link_source, best_link_p;
		std::vector<unsigned int> len_history, best_len_history;

		for (; p < str.size() && ms >= 0; p++, ms--)
		{
			char c = str[p];

			// explore known and matched nodes before erased nodes
			if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
			{
				t = m_Nodes[t].next.at(c);
				l++;
				len_history.push_back(m_Nodes[t].len);	// l <= m_Nodes[t].len. The matched string is a suffix of the longest string to this node.

				// putting this if in if makes sure that the end char is a match
				if (l > best_len)
				{
					best_node = t;
					best_len = l;
					best_pos = p;
					best_len_history = len_history;
				}

				best_link_depth = m_Nodes[m_Nodes[t].link].len;
				best_link_source = t;
				best_link_p = p;
			}
			else if (m_Nodes[t].next.find(ignored_char) != m_Nodes[t].next.end())
			{
				t = m_Nodes[t].next.at(ignored_char);
				l++;
				len_history.push_back(m_Nodes[t].len);	// l <= m_Nodes[t].len. The matched string is a suffix of the longest string to this node.

				// putting this if in if makes sure that the end char is a match
				if (l > best_len)
				{
					best_node = t;
					best_len = l;
					best_pos = p;
					best_len_history = len_history;
				}

				if (m_Nodes[m_Nodes[t].link].len >= best_link_depth)
				{
					best_link_depth = m_Nodes[m_Nodes[t].link].len;
					best_link_source = t;
					best_link_p = p;
				}
			}
			else	// fail to search
			{
				len_history.clear();
				if (best_link_depth > 0)
				{
					best_link_depth = 0;
					ms += p - best_link_p;
					p = best_link_p + 1;
					c = str[p];
					t = m_Nodes[best_link_source].link;
					l = m_Nodes[t].len;
				}
				unsigned int original_depth = m_EraseChainLen[t];
				while (t && m_Nodes[t].next.find(c) == m_Nodes[t].next.end()
					&& m_Nodes[t].next.find(ignored_char) == m_Nodes[t].next.end()
					&& l + ms > best_len)
				{
					// t = m_Nodes[t].link;
					t = m_DeepestLink[t];
					l = m_Nodes[t].len;
					p = p + m_EraseChainLen[t] - original_depth;
					ms = ms + original_depth - m_EraseChainLen[t];
					c = str[p];
					original_depth = m_EraseChainLen[t];
				}

				if (l + ms <= best_len)
				{
					break;
				}
			}
		}
		return best_len_history;		// Technically can autamatically apply std::move()
	}
	else
	{
		unsigned int t = start_node, l = start_len, ms = max_successor, p = start_pos;
		unsigned int best_link_depth = 0, best_link_source, best_link_p;
		std::vector<unsigned int> len_history;
		std::vector<unsigned int> len_history_dual;
		std::vector<char> c_history;
		std::vector<unsigned int> p_history;
		std::vector<std::vector<std::pair<unsigned int, unsigned int>>> host_count;
		std::vector<unsigned int> node_history;

		for (; p < str.size() && ms >= 0; p++, ms--)
		{
			char c = str[p];

			// explore known and matched nodes before erased nodes
			if (m_Nodes[t].next.find(c) != m_Nodes[t].next.end())
			{
				c_history.push_back(c);
				p_history.push_back(p);
				host_count.push_back(m_HostRangeCache[t]);
				node_history.push_back(t);

				best_link_depth = 0;
				t = m_Nodes[t].next.at(c);
				l++;
				len_history.push_back(m_Nodes[t].len);	// l <= m_Nodes[t].len. The matched string is a suffix of the longest string to this node.
				len_history_dual.push_back(m_Nodes[t].len);

				// putting this if in if makes sure that the end char is a match
				if (l > best_len)
				{
					best_node = t;
					best_len = l;
					best_pos = p;
				}
			}
			else if (m_Nodes[t].next.find(ignored_char) != m_Nodes[t].next.end())
			{
				c_history.push_back(ignored_char);
				p_history.push_back(p);
				host_count.push_back(m_HostRangeCache[t]);
				node_history.push_back(t);

				t = m_Nodes[t].next.at(ignored_char);
				l++;
				len_history.push_back(m_Nodes[t].len);	// l <= m_Nodes[t].len. The matched string is a suffix of the longest string to this node.
				len_history_dual.push_back(m_Nodes[t].len);

				// putting this if in if makes sure that the end char is a match
				if (l > best_len)
				{
					best_node = t;
					best_len = l;
					best_pos = p;
				}

				if (m_Nodes[m_Nodes[t].link].len > best_link_depth)
				{
					best_link_depth = m_Nodes[m_Nodes[t].link].len;
					best_link_source = t;
					best_link_p = p + 1;
				}
			}
			else	// fail to search
			{
				len_history.clear();
				if (best_link_depth > 0)
				{
					c = str[best_link_p];
					p = best_link_p;
					t = m_Nodes[best_link_source].link;
					l = m_Nodes[t].len;
				}

				while (t && m_Nodes[t].next.find(c) == m_Nodes[t].next.end()
					&& m_Nodes[t].next.find(ignored_char) == m_Nodes[t].next.end()
					&& l + ms > best_len)
				{
					t = m_Nodes[t].link;
					l = m_Nodes[t].len;
				}

				if (l + ms <= best_len)
				{
					break;
				}
			}
		}
		host_count.push_back(m_HostRangeCache[t]);
		node_history.push_back(t);

		if (p_history.empty())
			return len_history;

		std::cout << std::endl << str << std::endl;
		std::string spaces(p_history[0] - 1, ' ');
		std::cout << spaces << " ";
		for (auto& ph : p_history)
			std::cout << str[ph];
		std::cout << std::endl << spaces << " ";
		for (auto& c : c_history)
			std::cout << c;
		std::cout << std::endl << spaces;
		for (auto& th : node_history)
		{
			std::cout << m_Nodes[th].len;
		}
		std::cout << std::endl << spaces;
		/*for (auto& sc : host_count)
			std::cout << sc.coverlength();
		std::cout << std::endl;
		std::vector<std::vector<bool>> mat;
		for (auto& sc : host_count)
			mat.push_back(sc.get_individual_cover(0, 3));
		for (unsigned int i = 0; i < mat[0].size(); i++)
		{
			std::cout << spaces;
			for (unsigned int j = 0; j < mat.size(); j++)
			{
				std::cout << (mat[j][i] ? '*' : ' ');
			}
			std::cout << std::endl;
		}*/
		std::cout << "\nlen history" << std::endl;
		for (auto& i : len_history_dual)
		{
			std::cout << i << " ";
		}
		std::cout << std::endl;
		for (auto& i : node_history)
		{
			std::cout << i << " ";
		}
		std::cout << std::endl;
		return len_history;		// Technically can autamatically apply std::move()
	}
}

inline bool ChenAlignerCore::isNotCloned(const unsigned int& node) const
{
	return node < m_UnclonedLength;
}

SegmentTree::SegmentTree(size_t size)
{
	m_STMin.resize(4 * size, 0);
}

void print_segment(const std::string& str, const std::vector<std::pair<unsigned int, unsigned int>>& host_range, const unsigned int& l, const unsigned int& r)
{
	std::vector<bool> inlist;
	inlist.resize(4, false);
	for (auto& pr : host_range)
	{
		for (unsigned int i = pr.first; i <= pr.second; i++)
		{
			inlist[i] = true;
		}
	}

	std::string s1(l, ' ');
	std::string s2(r - l + 1, 'U');

	std::cout << "____" << std::endl;
	for (const auto& i : inlist)
	{
		if (i)
			std::cout << s1 << s2;
		std::cout << std::endl;
	}
	std::cout << "^^^^" << std::endl;
}

ChenAligner::ChenAligner(unsigned int alphabet_size) : m_Core(alphabet_size)
{
    // Nothing
}

ChenAligner::ChenAligner(unsigned int alphabet_size, char erase_marker) : m_Core(alphabet_size, erase_marker)
{
    // Nothing
}

void ChenAligner::GetWatermarksFromCache(const std::string& file_name)
{
    std::ifstream cache(file_name.c_str(), std::ios::binary);
    if (cache.is_open())
        std::cout << "File opened" << std::endl;
    else
        std::cout << "Fail to open file" << std::endl;
    m_Core.cache_acquireStructure(cache);

    init_finished = true;
}

void ChenAligner::DumpWatermarkStructureToCache(const std::string& file_name)
{
    if (!init_finished)
        throw std::runtime_error("No watermarks found. Please use GetWatermarks first!");

    std::ofstream cache(file_name.c_str(), std::ios::binary);
    m_Core.cache_storeStructure(cache);
}

void ChenAligner::MatchReads(ReadType read_type, std::string file_name, std::vector<unsigned int>& result, double sub_error_rate, bool use_parallel, unsigned int thread_count)
{
    if (!init_finished)
        throw std::runtime_error("No watermarks found. Please use GetWatermarks first!");

    if (use_parallel)
        MatchReadsParallel(read_type, file_name, result, sub_error_rate, thread_count);
    else
        MatchReadsNormal(read_type, file_name, result, sub_error_rate);
}

void ChenAligner::MatchReads(ReadType read_type, std::string file_name, std::string result_filename, double sub_error_rate, unsigned int thread_count)
{
    std::ifstream file(file_name.c_str());
    std::string read;
    unsigned int i = 0;
    if (!file.is_open())
    {
		std::cout << "reads file not opened" << std::endl;
        throw std::runtime_error("reads file not opened");
    }
    std::cout << "start matching" << std::endl;
    m_Core.parallel_Init(static_cast<unsigned int>(read.length() * sub_error_rate),
        static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
    switch (read_type)
    {
    case ReadType::FASTQ:
        while (std::getline(file, read))
        {
            std::getline(file, read);
            m_Core.parallel_addRead(read, static_cast<unsigned int>(read.length() * sub_error_rate),
                static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
            std::getline(file, read);
            std::getline(file, read);
            i++;
            if (i % 10000 == 0)
                std::cout << i << " reads set" << std::endl;
        }
        break;
    case ReadType::FASTA:
        while (std::getline(file, read))
        {
            std::getline(file, read);
            m_Core.parallel_addRead(read, static_cast<unsigned int>(read.length() * sub_error_rate),
                static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
            i++;
            if (i % 10000 == 0)
                std::cout << i << " reads set" << std::endl;
        }
        break;
    case ReadType::PURE:
        while (std::getline(file, read))
        {
            m_Core.parallel_addRead(read, static_cast<unsigned int>(read.length() * sub_error_rate),
                static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
            i++;
            if (i % 10000 == 0)
                std::cout << i << " reads set" << std::endl;
        }
        break;
    default:
        break;
    }
	auto start = std::chrono::high_resolution_clock::now();
    m_Core.parallel_terminateReadTransfer(thread_count);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Runtime: " << duration << std::endl;

    std::ofstream log_file(LOG_FILE_NAME, std::ios::app);
    log_file << result_filename << " TimeUse: " << duration << std::endl;
    log_file.close();

    std::cout << "Writing result to file" << std::endl;
    m_Core.parallel_getResults(result_filename, OUTFILE_FORMAT::ALL);
    std::cout << "Finish writing" << std::endl;
}

void ChenAligner::MatchReadsNormal(ReadType read_type, std::string file_name, std::vector<unsigned int>& result, double sub_error_rate)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::ifstream file(file_name.c_str());
    std::string read;
    unsigned int i = 0;
    if (!file.is_open())
    {
        throw std::runtime_error("reads file not opened");
    }
    std::cout << "start matching" << std::endl;
    switch (read_type)
    {
    case ReadType::FASTQ:
        while (std::getline(file, read))
        {
            std::getline(file, read);
            auto r = m_Core.alignRead(read, static_cast<unsigned int>(read.length() * sub_error_rate),
                static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
            if (r.empty())
            {
                result.push_back(0);
                std::cout << 0 << std::endl;
            }
            else
            {
                result.push_back(r[0].first);
                std::cout << r[0].first << std::endl;
            }
            std::getline(file, read);
            std::getline(file, read);
            i++;
            if (i % 1000 == 0)
                std::cout << i << " reads matched" << std::endl;
            //system("pause");
        }
        break;
    case ReadType::FASTA:
        while (std::getline(file, read))
        {
            std::getline(file, read);
            auto r = m_Core.alignRead(read, static_cast<unsigned int>(read.length() * sub_error_rate),
                static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
            if (r.empty())
                result.push_back(0);
            else
                result.push_back(r[0].first);
            i++;
            if (i % 1000 == 0)
                std::cout << i << " reads matched" << std::endl;
        }
        break;
    case ReadType::PURE:
        while (std::getline(file, read))
        {
            auto r = m_Core.alignRead(read, static_cast<unsigned int>(read.length() * sub_error_rate),
                static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
            if (r.empty())
                result.push_back(0);
            else
                result.push_back(r[0].first);
            i++;
            if (i % 1000 == 0)
                std::cout << i << " reads matched" << std::endl;
        }
        break;
    default:
        break;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Runtime: " << duration << std::endl;
}

void ChenAligner::MatchReadsParallel(ReadType read_type, std::string file_name, std::vector<unsigned int>& result, double sub_error_rate, unsigned int thread_count)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::ifstream file(file_name.c_str());
    std::string read;
    unsigned int i = 0;
    if (!file.is_open())
    {
        throw std::runtime_error("reads file not opened");
    }
    std::cout << "start matching" << std::endl;
    m_Core.parallel_Init(static_cast<unsigned int>(read.length() * sub_error_rate),
        static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
    switch (read_type)
    {
    case ReadType::FASTQ:
        while (std::getline(file, read))
        {
            std::getline(file, read);
            m_Core.parallel_addRead(read, static_cast<unsigned int>(read.length() * sub_error_rate),
                static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
            std::getline(file, read);
            std::getline(file, read);
            i++;
            if (i % 10000 == 0)
                std::cout << i << " reads set" << std::endl;
        }
        break;
    case ReadType::FASTA:
        while (std::getline(file, read))
        {
            std::getline(file, read);
            m_Core.parallel_addRead(read, static_cast<unsigned int>(read.length() * sub_error_rate),
                static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
            i++;
            if (i % 10000 == 0)
                std::cout << i << " reads set" << std::endl;
        }
        break;
    case ReadType::PURE:
        while (std::getline(file, read))
        {
            m_Core.parallel_addRead(read, static_cast<unsigned int>(read.length() * sub_error_rate),
                static_cast<unsigned int>(read.length() * INDEL_RATE), static_cast<unsigned int>(read.length() * LEN_DIFF_ACCEPTED));
            i++;
            if (i % 10000 == 0)
                std::cout << i << " reads set" << std::endl;
        }
        break;
    default:
        break;
    }
    m_Core.parallel_terminateReadTransfer(thread_count);
    m_Core.parallel_getResults(result);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Runtime: " << duration << std::endl;
}

void ChenAligner::GetWatermarksFromFile(const std::string& file_name)
{
    std::ifstream watermarks(file_name.c_str());
    if (!watermarks.is_open())
    {
        std::cout << "watermark file not opened" << std::endl;
        return;
    }

    std::cout << "start reading watermarks" << std::endl;
    std::string watermark;
    unsigned int i = 1;
    while (std::getline(watermarks, watermark))
    {
        m_Core.addString(watermark, i);
        i++;
    }
    watermarks.close();
    std::cout << "initializing" << std::endl;
    m_Core.build();
    std::cout << "finish initialization" << std::endl;

    init_finished = true;
}

std::string read2filename(const std::string& dna_sequence)
{
	if (dna_sequence.length() != 150)
		throw std::invalid_argument("DNA sequence must be exactly 150 characters long.");

	// Define the mapping for individual characters
	int char_to_num[128];  // Enough space for ASCII characters
	char_to_num['A'] = 0;
	char_to_num['C'] = 1;
	char_to_num['G'] = 2;
	char_to_num['T'] = 3;

	// Output character set (64 characters)
	std::string output_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";
	std::string result;

	// Process every three characters
	for (int i = 0; i < dna_sequence.length(); i += 3) {
		std::string chunk = dna_sequence.substr(i, 3);
		if (chunk.length() != 3)
			throw std::invalid_argument("DNA sequence chunks must divide evenly by 3.");

		// Calculate the index from three characters
		int index = 0;
		for (char ch : chunk) {
			index = index * 4 + char_to_num[ch];
		}

		// Map index to output character
		result.push_back(output_chars[index % 64]);
	}

	return result;
}