export interface Sponsor {
  id: string
  name: string
  logo: string
  website: string
  description: string
  tier: 'platinum' | 'gold' | 'silver'
}

// 临时数据 - 可以从API或其他数据源获取
export const mockSponsors: Sponsor[] = [
  {
    id: 'illumina',
    name: 'Illumina',
    logo: 'https://images.unsplash.com/photo-1581094794329-c8112a89af12?w=200&h=200&fit=crop&crop=faces',
    website: 'https://www.illumina.com',
    description: '全球领先的基因测序技术公司，为生物信息学研究提供创新的测序解决方案。',
    tier: 'platinum'
  },
  {
    id: 'thermo-fisher',
    name: 'Thermo Fisher Scientific',
    logo: 'https://images.unsplash.com/photo-1581094794329-c8112a89af12?w=200&h=200&fit=crop&crop=faces',
    website: 'https://www.thermofisher.com',
    description: '科学服务领域的领导者，提供分析仪器、实验室设备、试剂和耗材。',
    tier: 'platinum'
  },
  {
    id: 'qiagen',
    name: 'QIAGEN',
    logo: 'https://images.unsplash.com/photo-1581094794329-c8112a89af12?w=200&h=200&fit=crop&crop=faces',
    website: 'https://www.qiagen.com',
    description: '专注于样本制备和分子诊断技术的公司，为生命科学研究提供完整解决方案。',
    tier: 'gold'
  },
  {
    id: 'broad-institute',
    name: 'Broad Institute',
    logo: 'https://images.unsplash.com/photo-1581094794329-c8112a89af12?w=200&h=200&fit=crop&crop=faces',
    website: 'https://www.broadinstitute.org',
    description: '麻省理工学院和哈佛大学合作的生物医学和基因组学研究中心。',
    tier: 'gold'
  },
  {
    id: 'docker',
    name: 'Docker',
    logo: 'https://images.unsplash.com/photo-1581094794329-c8112a89af12?w=200&h=200&fit=crop&crop=faces',
    website: 'https://www.docker.com',
    description: '容器化技术的领导者，为生物信息学软件部署提供标准化解决方案。',
    tier: 'silver'
  },
  {
    id: 'microsoft',
    name: 'Microsoft',
    logo: 'https://images.unsplash.com/photo-1581094794329-c8112a89af12?w=200&h=200&fit=crop&crop=faces',
    website: 'https://www.microsoft.com',
    description: '通过Azure云平台为生物信息学研究提供强大的计算资源和AI服务。',
    tier: 'silver'
  }
]

export async function getAllSponsors(): Promise<Sponsor[]> {
  // 按赞助级别排序
  const tierOrder = { platinum: 3, gold: 2, silver: 1 }
  return mockSponsors.sort((a, b) => tierOrder[b.tier] - tierOrder[a.tier])
}

export function getSponsorsByTier(sponsors: Sponsor[], tier: Sponsor['tier']): Sponsor[] {
  return sponsors.filter(sponsor => sponsor.tier === tier)
}

export function getTierDisplayName(tier: Sponsor['tier']): string {
  switch (tier) {
    case 'platinum':
      return '白金赞助商'
    case 'gold':
      return '黄金赞助商'
    case 'silver':
      return '银牌赞助商'
    default:
      return '赞助商'
  }
}

export function getTierColor(tier: Sponsor['tier']): string {
  switch (tier) {
    case 'platinum':
      return 'from-gray-300 to-gray-500'
    case 'gold':
      return 'from-yellow-400 to-yellow-600'
    case 'silver':
      return 'from-gray-400 to-gray-600'
    default:
      return 'from-blue-400 to-blue-600'
  }
}